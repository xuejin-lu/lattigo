package ckks

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

func experimentalTrustTestParameters(t *testing.T) Parameters {
	t.Helper()
	// CNN's L26 chain shape and L15 entry, at small N for a zero-Bootstrap test.
	logQ := []int{42, 32, 32, 32}
	for i := 0; i < 12; i++ {
		logQ = append(logQ, 32)
	}
	for i := 0; i < 11; i++ {
		logQ = append(logQ, 42)
	}
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: 10, LogQ: logQ, LogP: []int{42, 42, 42, 42, 42}, LogDefaultScale: 32,
	})
	require.NoError(t, err)
	return params
}

func TestFastQPExperimentalTrustL15AcrossEvaluators(t *testing.T) {
	params := experimentalTrustTestParameters(t)
	const level, logSlots, baby, giant = 15, 8, 1, 64
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	gkA := kgen.GenGaloisKeyNew(params.GaloisElement(baby), sk)
	gkB := kgen.GenGaloisKeyNew(params.GaloisElement(giant), sk)
	a := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gkA))
	b := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gkB))
	require.NotSame(t, a.fastQPCompat, b.fastQPCompat)
	require.True(t, a.usesFastQPExperimentalTrust())
	require.True(t, b.usesFastQPExperimentalTrust())
	encoder := NewEncoder(params)
	values := make([]complex128, 1<<logSlots)
	for i := range values {
		values[i] = complex(float64(i%17-8)/32, float64(i%7-3)/64)
	}
	pt := NewPlaintext(params, level)
	pt.LogDimensions = ring.Dimensions{Cols: logSlots}
	require.NoError(t, encoder.Encode(values, pt))
	input, err := rlwe.NewEncryptor(params, pk).EncryptNew(pt)
	require.NoError(t, err)
	inputBefore := input.CopyNew()
	levelP := params.MaxLevelP()
	ringQP := params.RingQP().AtLevel(level, levelP)
	decomposition := make([]ringqp.Poly, params.BaseRNSDecompositionVectorSize(level, levelP))
	for i := range decomposition {
		decomposition[i] = ringQP.NewPoly()
	}
	a.DecomposeNTT(level, levelP, levelP+1, input.Value[1], input.IsNTT, decomposition)
	rotatedQP := rlwe.NewElementExtended(params, 1, level, levelP)
	require.NoError(t, a.AutomorphismHoistedLazy(level, input, decomposition, params.GaloisElement(baby), rotatedQP))

	// The same actual full-QP Montgomery mask multiply/add used by CNN512.
	mask := ringQP.NewPoly()
	maskMeta := *pt.MetaData
	maskMeta.Scale = rlwe.NewScale(params.Q()[level])
	maskMeta.IsMontgomery = true
	maskValues := make([]float64, len(values))
	for i := range maskValues {
		maskValues[i] = 0.25
	}
	require.NoError(t, encoder.Embed(maskValues, &maskMeta, mask))
	acc := rlwe.NewElementExtended(params, 1, level, levelP)
	for term := 0; term < 2; term++ {
		for component := range acc.Value {
			ringQP.MulCoeffsMontgomeryThenAdd(mask, rotatedQP.Value[component], acc.Value[component])
		}
	}
	accBefore := acc.CopyNew()
	strict := NewEvaluator(params, nil)
	strict.fastQPCompat.experimentalTrust = false
	strictOut := params.RingQ().AtLevel(level).NewPoly()
	strict.ModDownQPtoQNTT(level, levelP, acc.Value[1].Q, acc.Value[1].P, strictOut)
	require.ErrorContains(t, strict.fastQPCompat.peek(), "high-level dormant Q rows are not authorized")
	require.False(t, b.fastQPCompat.isArmed())
	b.ModDownQPtoQNTT(level, levelP, acc.Value[1].Q, acc.Value[1].P, acc.Value[1].Q)
	require.NoError(t, b.fastQPCompat.peek())
	ksOut := rlwe.NewElementExtended(params, 1, level, levelP)
	require.NoError(t, b.GadgetProductLazy(level, acc.Value[1].Q, &gkB.GadgetCiphertext, ksOut))
	ringQP.Add(ksOut.Value[0], acc.Value[0], ksOut.Value[0])
	final := rlwe.NewElementExtended(params, 1, level, levelP)
	for component := range final.Value {
		ringQP.AutomorphismNTTWithIndex(ksOut.Value[component], b.AutomorphismIndex(params.GaloisElement(giant)), final.Value[component])
	}
	output := NewCiphertext(params, 1, level)
	*output.MetaData = *pt.MetaData
	output.Scale = pt.Scale.Mul(maskMeta.Scale)
	for component := range output.Value {
		b.ModDownQPtoQNTT(level, levelP, final.Value[component].Q, final.Value[component].P, output.Value[component])
	}
	require.NoError(t, b.fastQPCompat.peek())
	require.True(t, fastQPolyIsZero(output.Value[1], level+1))
	require.False(t, b.fastQPCompat.isArmed(), "trust must not fake an armed history")
	require.Empty(t, b.fastQPCompat.trusted, "trust must not hash or certify ciphertexts")
	require.Equal(t, uint64(1), a.fastQPCompat.counterSnapshot().keySwitchElisions)
	require.Equal(t, uint64(1), b.fastQPCompat.counterSnapshot().gadgetProductElisions)
	require.True(t, input.Equal(inputBefore))
	require.Equal(t, accBefore.Value[0], acc.Value[0], "ModDown of c1 must not overwrite c0")

	want := fastRotateTestOracle(values, baby+giant)
	for i := range want {
		want[i] *= 0.5 // two mask terms, each 0.25
	}
	decryptor := rlwe.NewDecryptor(params, sk)
	checkDecoded := func(label string, ct *rlwe.Ciphertext) {
		t.Helper()
		got := make([]complex128, len(want))
		require.NoError(t, encoder.Decode(decryptor.DecryptNew(ct), got))
		var maxError float64
		for i := range want {
			maxError = math.Max(maxError, cmplx.Abs(got[i]-want[i]))
		}
		t.Logf("%s: max complex error %.12g (limit 1e-6)", label, maxError)
		require.LessOrEqual(t, maxError, 1e-6)
	}
	checkDecoded("L15 actual accumulated ModDown", output)

	// A third, fresh evaluator must route complete-Q Rescale correctly too.
	c := NewEvaluator(params, nil)
	before := output.CopyNew()
	oracle := NewCiphertext(params, 1, level-1)
	buffer := params.RingQ().AtLevel(level).NewPoly()
	for component := range oracle.Value {
		params.RingQ().AtLevel(level).DivRoundByLastModulusManyNTT(1, before.Value[component], buffer, oracle.Value[component])
	}
	require.NoError(t, c.Rescale(output, output))
	require.Equal(t, oracle.Value, output.Value, "complete-Q Rescale must equal the full-RNS oracle on every row")
	require.Equal(t, uint64(1), c.fastQPCompat.counterSnapshot().fullQRescales)
	require.Equal(t, level-1, output.Level())
	require.True(t, before.Scale.Div(rlwe.NewScale(params.Q()[level])).Equal(output.Scale))
	checkDecoded("L14 actual Rescale", output)

	// Complete-Q trust must not steal the compact Q-prefix route.
	compact := before.CopyNew()
	for component := range compact.Value {
		for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
			compact.Value[component].Coeffs[row] = nil
		}
	}
	d := NewEvaluator(params, nil)
	require.NoError(t, d.Rescale(compact, compact))
	require.Zero(t, d.fastQPCompat.counterSnapshot().fullQRescales)
	for row := fastcore.MaxQPrefixWidth; row <= compact.Level(); row++ {
		require.Nil(t, compact.Value[0].Coeffs[row])
	}
	checkDecoded("L14 shared Q-prefix Rescale", compact)
}

func TestFastQPExperimentalTrustRetainsStructuralAndZeroGuards(t *testing.T) {
	params := experimentalTrustTestParameters(t)
	const level = 15
	levelP := params.MaxLevelP()
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	gk := kgen.GenGaloisKeyNew(params.GaloisElement(64), sk)
	for _, name := range []string{"missing-Q", "missing-P", "missing-output", "nonzero-P", "bad-level"} {
		t.Run(name, func(t *testing.T) {
			eval := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
			qp := rlwe.NewElementExtended(params, 1, level, levelP)
			out := params.RingQ().AtLevel(level).NewPoly()
			requestedLevel := level
			switch name {
			case "missing-Q":
				qp.Value[0].Q.Coeffs[4] = nil
			case "missing-P":
				qp.Value[0].P.Coeffs[1] = nil
			case "missing-output":
				out.Coeffs[4] = nil
			case "nonzero-P":
				qp.Value[0].P.Coeffs[0][0] = 1
			case "bad-level":
				requestedLevel = params.MaxLevel() + 1
			}
			// Snapshot deliberately malformed storage without passing it through
			// Poly.CopyNew's new compact-shape validation. Oracle and guards stay unchanged.
			before := &ring.Poly{Coeffs: out.Coeffs.CopyNew()}
			// Matrix.CopyNew normalizes nil rows to empty slices; retain their shape.
			for row := range out.Coeffs {
				if out.Coeffs[row] == nil {
					before.Coeffs[row] = nil
				}
			}
			eval.ModDownQPtoQNTT(requestedLevel, levelP, qp.Value[0].Q, qp.Value[0].P, out)
			require.Error(t, eval.fastQPCompat.peek())
			require.Equal(t, *before, out, "failed ModDown must leave output untouched")
			err := eval.GadgetProductLazy(level, qp.Value[1].Q, &gk.GadgetCiphertext, qp)
			require.ErrorContains(t, err, "previous compatibility operation failed")
		})
	}
	eval := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
	qp := rlwe.NewElementExtended(params, 1, level, levelP)
	qp.Value[1].Q.Coeffs[0][0] = 1
	require.ErrorContains(t, eval.GadgetProductLazy(level, qp.Value[1].Q, &gk.GadgetCiphertext, qp), "source polynomial is nonzero")
}

func TestFastQPExperimentalTrustSkipsOnlyNumericalConsistency(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	input := NewCiphertext(params, 1, 5)
	input.Value[0].Coeffs[4][0] = 1 // deliberately inconsistent with the zero prefix
	trust := NewEvaluator(params, nil)
	strict := NewEvaluator(params, nil)
	strict.fastQPCompat.experimentalTrust = false
	require.NoError(t, trust.verifyFastFullQCiphertext(input))
	require.ErrorContains(t, strict.verifyFastFullQCiphertext(input), "inconsistent")
	input.Value[0].Coeffs[4] = nil
	require.Error(t, trust.verifyFastFullQCiphertext(input), "trust must not invent missing rows")
	require.False(t, trust.fastQPCompat.isArmed())
}
