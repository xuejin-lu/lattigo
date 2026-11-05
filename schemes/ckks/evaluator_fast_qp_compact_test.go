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

type compactQPRescaleRecorder struct {
	fastcore.RescaleCore
	rows     []int
	divisors []uint64
}

func (r *compactQPRescaleRecorder) ApplyRows(q *ring.Ring, in, out *rlwe.Ciphertext, count, rows int) error {
	r.rows = append(r.rows, rows)
	r.divisors = append(r.divisors, q.SubRings[in.Level()].Modulus)
	return r.RescaleCore.ApplyRows(q, in, out, count, rows)
}

func TestFastQPCompactAcrossEvaluatorsUsesSharedCores(t *testing.T) {
	params := experimentalTrustTestParameters(t)
	const level, logSlots = 15, 8
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	gk := kgen.GenGaloisKeyNew(params.GaloisElement(1), sk)
	encoder := NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	values := make([]complex128, 1<<logSlots)
	for i := range values {
		values[i] = complex(float64(i%17-8)/32, float64(i%7-3)/64)
	}
	pt := NewPlaintext(params, level)
	pt.LogDimensions = ring.Dimensions{Cols: logSlots}
	require.NoError(t, encoder.Encode(values, pt))
	native, err := rlwe.NewEncryptor(params, kgen.GenPublicKeyNew(sk)).EncryptNew(pt)
	require.NoError(t, err)
	input, err := NewEvaluator(params, nil).AddNew(native, native)
	require.NoError(t, err)
	before := input.CopyNew()
	levelP := params.MaxLevelP()
	rqp := params.RingQP().AtLevel(level, levelP)
	for _, strict := range []bool{false, true} {
		for _, alias := range []bool{false, true} {
			t.Run(map[bool]string{true: "strict", false: "trust"}[strict]+map[bool]string{true: "/alias", false: "/separate"}[alias], func(t *testing.T) {
				a := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
				a.fastQPCompat.experimentalTrust = !strict
				recorder := &recordingPublicAutomorphismCore{AutomorphismCore: a.automorphismCore}
				a.automorphismCore = recorder
				decomp := make([]ringqp.Poly, params.BaseRNSDecompositionVectorSize(level, levelP))
				for i := range decomp {
					decomp[i] = rqp.NewPoly()
				}
				a.DecomposeNTT(level, levelP, levelP+1, input.Value[1], true, decomp)
				qp := rlwe.NewElementExtended(params, 1, level, levelP)
				if alias {
					qp.Value[0].Q = input.Value[0]
				}
				require.NoError(t, a.AutomorphismHoistedLazy(level, input, decomp, params.GaloisElement(1), qp))
				require.Equal(t, []int{4}, recorder.rows)
				require.Equal(t, uint64(1), a.fastQPCompat.counterSnapshot().compactHoistedAutomorphisms)
				require.Equal(t, uint64(1), a.fastQPCompat.counterSnapshot().keySwitchElisions)
				require.True(t, input.Equal(before))
				for _, component := range qp.Value {
					for row := 4; row <= level; row++ {
						require.Empty(t, component.Q.Coeffs[row])
					}
					require.True(t, fastQPolyIsZero(component.P, levelP+1))
				}
				mask := rqp.NewPoly()
				meta := *pt.MetaData
				meta.Scale = rlwe.NewScale(params.Q()[level])
				meta.IsMontgomery = true
				maskValues := make([]float64, len(values))
				for i := range maskValues {
					maskValues[i] = 0.25
				}
				require.NoError(t, encoder.Embed(maskValues, &meta, mask))
				acc := rlwe.NewElementExtended(params, 1, level, levelP)
				for component := range acc.Value {
					rqp.MulCoeffsMontgomeryThenAdd(mask, qp.Value[component], acc.Value[component])
				}
				b := NewEvaluator(params, nil)
				b.fastQPCompat.experimentalTrust = !strict
				out := NewCiphertext(params, 1, level)
				*out.MetaData = *input.MetaData
				out.Scale = input.Scale.Mul(meta.Scale)
				if alias {
					for component := range out.Value {
						out.Value[component] = acc.Value[component].Q
					}
				}
				for component := range out.Value {
					b.ModDownQPtoQNTT(level, levelP, acc.Value[component].Q, acc.Value[component].P, out.Value[component])
				}
				require.NoError(t, b.fastQPCompat.peek())
				require.False(t, b.fastQPCompat.isArmed())
				require.Equal(t, uint64(2), b.fastQPCompat.counterSnapshot().compactModDowns)
				check := func(ct *rlwe.Ciphertext) {
					t.Helper()
					got := make([]complex128, len(values))
					require.NoError(t, encoder.Decode(decryptor.DecryptNew(ct), got))
					for i, v := range got {
						require.False(t, math.IsNaN(real(v)) || math.IsNaN(imag(v)) || math.IsInf(cmplx.Abs(v), 0))
						require.LessOrEqual(t, cmplx.Abs(v-0.5*values[(i+1)%len(values)]), 1e-6)
					}
					for _, p := range ct.Value {
						for row := 4; row <= ct.Level(); row++ {
							require.Empty(t, p.Coeffs[row])
						}
					}
				}
				check(out)
				c := NewEvaluator(params, nil)
				c.fastQPCompat.experimentalTrust = !strict
				rescaleRecorder := &compactQPRescaleRecorder{RescaleCore: c.rescaleCore}
				c.rescaleCore = rescaleRecorder
				target := NewCiphertext(params, 1, level-1)
				if alias {
					target = out
				}
				scale := out.Scale.Div(rlwe.NewScale(params.Q()[level]))
				require.NoError(t, c.Rescale(out, target))
				require.Equal(t, level-1, target.Level())
				require.True(t, target.Scale.Equal(scale))
				require.Equal(t, []int{4}, rescaleRecorder.rows)
				require.Equal(t, []uint64{params.Q()[level]}, rescaleRecorder.divisors)
				require.Zero(t, c.fastQPCompat.counterSnapshot().fullQRescales)
				check(target)
			})
		}
	}
}

func TestFastQPCompactRejectsMissingRowsAndDefersErrors(t *testing.T) {
	params := experimentalTrustTestParameters(t)
	const level = 15
	for _, kind := range []string{"zero-width", "too-narrow", "hole", "nonzero-P", "missing-P", "output-hole"} {
		t.Run(kind, func(t *testing.T) {
			eval := NewEvaluator(params, nil)
			q := fastcore.NewCompactCiphertext(params, 1, level).Value[0]
			p := params.RingP().NewPoly()
			out := params.RingQ().AtLevel(level).NewPoly()
			switch kind {
			case "zero-width":
				for row := range q.Coeffs {
					q.Coeffs[row] = nil
				}
			case "too-narrow":
				q.Coeffs[3] = nil
			case "hole":
				q.Coeffs[4] = nil
				q.Coeffs[5] = make([]uint64, params.N())
			case "nonzero-P":
				p.Coeffs[0][0] = 1
			case "missing-P":
				p.Coeffs[1] = nil
			case "output-hole":
				out.Coeffs[4] = nil
			}
			before := out.Coeffs.CopyNew()
			for row := range before {
				if out.Coeffs[row] == nil {
					before[row] = nil
				}
			}
			eval.ModDownQPtoQNTT(level, params.MaxLevelP(), q, p, out)
			require.Error(t, eval.fastQPCompat.peek())
			require.Equal(t, before, out.Coeffs)
			require.Zero(t, eval.fastQPCompat.counterSnapshot().compactModDowns)
			input := fastcore.NewCompactCiphertext(params, 1, level)
			target := NewCiphertext(params, 1, level-1)
			targetBefore := target.CopyNew()
			require.ErrorContains(t, eval.Rescale(input, target), "previous compatibility operation failed")
			require.True(t, target.Equal(targetBefore))
			require.NoError(t, eval.fastQPCompat.peek())
		})
	}
}

func TestFastQPCompactHoistedGuardsAndFullBufferCompatibility(t *testing.T) {
	params := experimentalTrustTestParameters(t)
	const level = 15
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	galEl := params.GaloisElement(1)
	gk := kgen.GenGaloisKeyNew(galEl, sk)
	levelP := params.MaxLevelP()
	for _, kind := range []string{"zero-width", "too-narrow", "nonzero-c1", "nonzero-decomposition", "output-hole"} {
		t.Run(kind, func(t *testing.T) {
			eval := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
			recorder := &recordingPublicAutomorphismCore{AutomorphismCore: eval.automorphismCore}
			eval.automorphismCore = recorder
			in := fastcore.NewCompactCiphertext(params, 1, level)
			decomp := make([]ringqp.Poly, params.BaseRNSDecompositionVectorSize(level, levelP))
			for i := range decomp {
				decomp[i] = params.RingQP().AtLevel(level, levelP).NewPoly()
			}
			out := rlwe.NewElementExtended(params, 1, level, levelP)
			before := out.CopyNew()
			switch kind {
			case "zero-width":
				for row := range in.Value[0].Coeffs {
					in.Value[0].Coeffs[row] = nil
				}
			case "too-narrow":
				in.Value[0].Coeffs[3] = nil
			case "nonzero-c1":
				in.Value[1].Coeffs[1][0] = 1
			case "nonzero-decomposition":
				decomp[0].P.Coeffs[0][0] = 1
			case "output-hole":
				out.Value[0].Q.Coeffs[4] = nil
				before.Value[0].Q.Coeffs[4] = nil
			}
			require.Error(t, eval.AutomorphismHoistedLazy(level, in, decomp, galEl, out))
			require.True(t, out.Equal(before))
			require.Empty(t, recorder.rows)
		})
	}
	// Preserve the old full-Q lane's observable writes through a retained view.
	eval := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
	in := NewCiphertext(params, 1, level)
	in.Value[0].Coeffs[0][0] = 7
	decomp := make([]ringqp.Poly, params.BaseRNSDecompositionVectorSize(level, levelP))
	for i := range decomp {
		decomp[i] = params.RingQP().AtLevel(level, levelP).NewPoly()
	}
	out := rlwe.NewElementExtended(params, 1, level, levelP)
	view := out.Value[0].Q
	require.NoError(t, eval.AutomorphismHoistedLazy(level, in, decomp, galEl, out))
	require.True(t, view.Equal(&out.Value[0].Q))
	require.False(t, fastQPolyIsZero(view, level+1))
	require.Zero(t, eval.fastQPCompat.counterSnapshot().compactHoistedAutomorphisms)
}
