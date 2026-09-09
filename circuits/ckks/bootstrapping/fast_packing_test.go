package bootstrapping

import (
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func fastPackingParameters(t *testing.T, factorTwo bool) (Parameters, ckks.Parameters, ckks.Parameters) {
	t.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(50, 64)
	moduli, err := generator.NextAlternatingPrimes(5)
	require.NoError(t, err)
	n1, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 4, Q: moduli, LogDefaultScale: 30})
	require.NoError(t, err)
	n2 := n1
	if factorTwo {
		n2, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 5, Q: moduli, LogDefaultScale: 30})
		require.NoError(t, err)
	}
	maxLogSlots := 3
	if factorTwo {
		maxLogSlots = 4
	}
	params := Parameters{
		ResidualParameters:      n1,
		BootstrappingParameters: n2,
		SlotsToCoeffsParameters: dft.MatrixLiteral{LogSlots: maxLogSlots},
	}
	return params, n1, n2
}

func newFastPackingCiphertext(params ckks.Parameters, level, logSlots int, seed uint64) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, level)
	ct.IsNTT = true
	ct.IsMontgomery = true
	ct.IsBatched = true
	ct.IsBitReversed = true
	ct.Scale = rlwe.NewScale(1 << 30)
	ct.LogDimensions = ring.Dimensions{Rows: 0, Cols: logSlots}
	for d := 0; d <= 1; d++ {
		for limb, subring := range params.RingQ().SubRings[:2] {
			for i := range ct.Value[d].Coeffs[limb] {
				value := (seed + uint64(17*d+31*limb+i)) % subring.Modulus
				ct.Value[d].Coeffs[limb][i] = ring.MForm(value, subring.Modulus, subring.BRedConstant)
			}
		}
		for limb, subring := range params.RingQ().SubRings[2 : level+1] {
			for i := range ct.Value[d].Coeffs[limb+2] {
				ct.Value[d].Coeffs[limb+2][i] = ^uint64(0) - seed - uint64(d+limb+i)
				ct.Value[d].Coeffs[limb+2][i] %= subring.Modulus
			}
		}
	}
	return ct
}

func cloneFastPackingInputs(cts []rlwe.Ciphertext) []rlwe.Ciphertext {
	clones := make([]rlwe.Ciphertext, len(cts))
	for i := range cts {
		clones[i] = *cts[i].CopyNew()
	}
	return clones
}

func requireFastPackingQ01Equal(t *testing.T, want, got []rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, len(want), len(got))
	for i := range want {
		require.Equal(t, want[i].Degree(), got[i].Degree(), "ciphertext %d degree", i)
		require.Equal(t, *want[i].MetaData, *got[i].MetaData, "ciphertext %d metadata", i)
		for d := 0; d <= 1; d++ {
			require.Equal(t, want[i].Value[d].Coeffs[:2], got[i].Value[d].Coeffs[:2], "ciphertext %d component %d", i, d)
		}
	}
}

func standardPackReference(t *testing.T, params Parameters, cts []rlwe.Ciphertext, ctxt packingContext) []rlwe.Ciphertext {
	t.Helper()
	standard := Evaluator{Parameters: params}
	level := cts[0].Level()
	var xPow []ring.Poly
	if ctxt.Params.N() == params.ResidualParameters.N() && params.ResidualParameters.N() != params.BootstrappingParameters.N() {
		xPow = rlwe.GenXPow2NTT(ctxt.Params.RingQ().AtLevel(level), params.BootstrappingParameters.LogN(), false)
	} else {
		xPow = rlwe.GenXPow2NTT(ctxt.Params.RingQ().AtLevel(level), ctxt.Params.LogN(), false)
	}
	got, err := standard.pack(cloneFastPackingInputs(cts), ctxt, xPow)
	require.NoError(t, err)
	return got
}

func standardUnpackReference(t *testing.T, params Parameters, ct *rlwe.Ciphertext, ctxt packingContext, inverse bool) []rlwe.Ciphertext {
	t.Helper()
	standard := Evaluator{Parameters: params}
	level := ct.Level()
	logN := ctxt.Params.LogN()
	if !inverse && ctxt.Params.N() == params.ResidualParameters.N() && params.ResidualParameters.N() != params.BootstrappingParameters.N() {
		logN = params.BootstrappingParameters.LogN()
	}
	xPow := rlwe.GenXPow2NTT(ctxt.Params.RingQ().AtLevel(level), logN, true)
	got, err := standard.unpack(ct.CopyNew(), ctxt, xPow)
	require.NoError(t, err)
	return got
}

func TestFastPackMatchesStandardForCounts(t *testing.T) {
	params, n1, _ := fastPackingParameters(t, false)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	for _, count := range []int{1, 2, 3, 5} {
		t.Run(string(rune('0'+count)), func(t *testing.T) {
			inputs := make([]rlwe.Ciphertext, count)
			for i := range inputs {
				inputs[i] = *newFastPackingCiphertext(n1, 4, 1, uint64(100+i*19))
			}
			ctx := &packingContext{Params: &params.BootstrappingParameters, LogMaxDimensions: params.LogMaxDimensions(), LogSlots: 1, NbPackedCTs: count}
			want := standardPackReference(t, params, inputs, *ctx)
			got, _, _, err := eval.PackAndSwitchN1ToN2(inputs)
			require.NoError(t, err)
			requireFastPackingQ01Equal(t, want, got)
		})
	}
}

func TestFastUnpackMatchesStandardAndRoundTrips(t *testing.T) {
	params, n1, _ := fastPackingParameters(t, false)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	inputs := make([]rlwe.Ciphertext, 3)
	for i := range inputs {
		inputs[i] = *newFastPackingCiphertext(n1, 4, 1, uint64(200+i*23))
	}
	packed, _, packedCtx, err := eval.PackAndSwitchN1ToN2(inputs)
	require.NoError(t, err)
	want := standardUnpackReference(t, params, &packed[0], *packedCtx, true)
	for i := range want {
		want[i].LogDimensions.Cols = inputs[0].LogSlots()
	}
	ctxCopy := *packedCtx
	got, err := eval.UnpackAndSwitchN2ToN1(packed, nil, &ctxCopy)
	require.NoError(t, err)
	requireFastPackingQ01Equal(t, want, got)
	for i := range got {
		require.Equal(t, inputs[i].LogSlots(), got[i].LogSlots())
	}

	packed, _, packedCtx, err = eval.PackAndSwitchN1ToN2(inputs)
	require.NoError(t, err)
	ctxCopy = *packedCtx
	roundTrip, err := eval.UnpackAndSwitchN2ToN1(packed, nil, &ctxCopy)
	require.NoError(t, err)
	wantRoundTrip := standardUnpackReference(t, params, &packed[0], *packedCtx, true)
	for i := range wantRoundTrip {
		wantRoundTrip[i].LogDimensions.Cols = inputs[0].LogSlots()
	}
	requireFastPackingQ01Equal(t, wantRoundTrip, roundTrip)
}

func TestFastPackingIgnoresDormantResiduesAndOwnsOutputs(t *testing.T) {
	params, n1, _ := fastPackingParameters(t, false)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	inputsA := []rlwe.Ciphertext{*newFastPackingCiphertext(n1, 4, 1, 301), *newFastPackingCiphertext(n1, 4, 1, 307), *newFastPackingCiphertext(n1, 4, 1, 313)}
	inputsB := cloneFastPackingInputs(inputsA)
	for d := range inputsB[1].Value {
		for limb := 2; limb <= inputsB[1].Level(); limb++ {
			for i := range inputsB[1].Value[d].Coeffs[limb] {
				inputsB[1].Value[d].Coeffs[limb][i] ^= uint64(0x5a5a + d + limb + i)
			}
		}
	}
	packedA, _, ctxA, err := eval.PackAndSwitchN1ToN2(inputsA)
	require.NoError(t, err)
	packedB, _, _, err := eval.PackAndSwitchN1ToN2(inputsB)
	require.NoError(t, err)
	requireFastPackingQ01Equal(t, packedA, packedB)

	first := packedA[0].Value[0].Coeffs[0][0]
	packedAgain, _, _, err := eval.PackAndSwitchN1ToN2(inputsA)
	require.NoError(t, err)
	packedA[0].Value[0].Coeffs[0][0] ^= 1
	require.Equal(t, first, packedAgain[0].Value[0].Coeffs[0][0])

	ctxCopy := *ctxA
	unpacked, err := eval.UnpackAndSwitchN2ToN1(packedAgain, nil, &ctxCopy)
	require.NoError(t, err)
	unpacked[0].Value[0].Coeffs[0][0] ^= 1
	require.NotEqual(t, unpacked[0].Value[0].Coeffs[0][0], unpacked[1].Value[0].Coeffs[0][0])
}

func TestFastPackingFactorTwoRingBoundary(t *testing.T) {
	params, n1, n2 := fastPackingParameters(t, true)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	inputs := make([]rlwe.Ciphertext, 5)
	for i := range inputs {
		inputs[i] = *newFastPackingCiphertext(n1, 4, 1, uint64(401+i*29))
	}
	got, ctxN1, ctxN2, err := eval.PackAndSwitchN1ToN2(inputs)
	require.NoError(t, err)
	require.NotNil(t, ctxN1)
	require.NotNil(t, ctxN2)
	require.Equal(t, n2.N(), got[0].N())

	packN1 := &packingContext{Params: &n1, LogMaxDimensions: n1.LogMaxDimensions(), LogSlots: 1, NbPackedCTs: len(inputs)}
	wantN1, err := (Evaluator{Parameters: params}).pack(cloneFastPackingInputs(inputs), *packN1, rlwe.GenXPow2NTT(n1.RingQ().AtLevel(4), n2.LogN(), false))
	require.NoError(t, err)
	wantN2Inputs := make([]rlwe.Ciphertext, len(wantN1))
	for i := range wantN1 {
		out := ckks.NewCiphertext(n2, 1, wantN1[i].Level())
		out.IsNTT, out.IsMontgomery = wantN1[i].IsNTT, wantN1[i].IsMontgomery
		rlwe.SwitchCiphertextRingDegreeNTT(wantN1[i].El(), nil, out.El())
		wantN2Inputs[i] = *out
	}
	packN2 := &packingContext{Params: &n2, LogMaxDimensions: params.LogMaxDimensions(), LogSlots: wantN2Inputs[0].LogSlots(), NbPackedCTs: len(wantN2Inputs)}
	want, err := (Evaluator{Parameters: params}).pack(wantN2Inputs, *packN2, rlwe.GenXPow2NTT(n2.RingQ().AtLevel(4), n2.LogN(), false))
	require.NoError(t, err)
	requireFastPackingQ01Equal(t, want, got)

	ctxN1Copy, ctxN2Copy := *ctxN1, *ctxN2
	gotUnpacked, err := eval.UnpackAndSwitchN2ToN1(got, &ctxN1Copy, &ctxN2Copy)
	require.NoError(t, err)
	require.Len(t, gotUnpacked, len(inputs))
	refN2Ctx, refN1Ctx := *packN2, *packN1
	var refN2Unpacked []rlwe.Ciphertext
	for i := range want {
		unpacked := standardUnpackReference(t, params, &want[i], refN2Ctx, true)
		refN2Unpacked = append(refN2Unpacked, unpacked...)
		refN2Ctx.NbPackedCTs -= len(unpacked)
	}
	var refUnpacked []rlwe.Ciphertext
	for i := range refN2Unpacked {
		out := ckks.NewCiphertext(n1, 1, refN2Unpacked[i].Level())
		out.IsNTT, out.IsMontgomery = refN2Unpacked[i].IsNTT, refN2Unpacked[i].IsMontgomery
		rlwe.SwitchCiphertextRingDegreeNTT(refN2Unpacked[i].El(), n2.RingQ(), out.El())
		unpacked := standardUnpackReference(t, params, out, refN1Ctx, true)
		refUnpacked = append(refUnpacked, unpacked...)
		refN1Ctx.NbPackedCTs -= len(unpacked)
	}
	for i := range refUnpacked {
		refUnpacked[i].LogDimensions.Cols = inputs[0].LogSlots()
	}
	requireFastPackingQ01Equal(t, refUnpacked, gotUnpacked)
	for i := range gotUnpacked {
		require.Equal(t, n1.N(), gotUnpacked[i].N())
		require.Equal(t, inputs[i].LogSlots(), gotUnpacked[i].LogSlots())
	}
}

func TestFastPackingValidation(t *testing.T) {
	params, n1, _ := fastPackingParameters(t, false)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	_, _, _, err = eval.PackAndSwitchN1ToN2(nil)
	require.Error(t, err)
	bad := *newFastPackingCiphertext(n1, 4, 1, 503)
	bad.IsNTT = false
	_, _, _, err = eval.PackAndSwitchN1ToN2([]rlwe.Ciphertext{bad})
	require.Error(t, err)
	bad = *newFastPackingCiphertext(n1, 4, 1, 509)
	bad.Resize(2, bad.Level())
	_, _, _, err = eval.PackAndSwitchN1ToN2([]rlwe.Ciphertext{bad})
	require.Error(t, err)

	ci, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 4, Q: n1.Q(), RingType: ring.ConjugateInvariant, LogDefaultScale: 30})
	require.NoError(t, err)
	ciParams := params
	ciParams.BootstrappingParameters = ci
	_, err = NewFastEvaluator(ciParams)
	require.Error(t, err)

	generator := ring.NewNTTFriendlyPrimesGenerator(50, 128)
	moduli, err := generator.NextAlternatingPrimes(5)
	require.NoError(t, err)
	small, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 4, Q: moduli, LogDefaultScale: 30})
	require.NoError(t, err)
	large, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 6, Q: moduli, LogDefaultScale: 30})
	require.NoError(t, err)
	tooWide := params
	tooWide.ResidualParameters = small
	tooWide.BootstrappingParameters = large
	tooWideEval, err := NewFastEvaluator(tooWide)
	require.NoError(t, err)
	_, _, _, err = tooWideEval.PackAndSwitchN1ToN2([]rlwe.Ciphertext{*newFastPackingCiphertext(small, 4, 1, 521)})
	require.Error(t, err)

	tooSmall := params
	tooSmall.ResidualParameters = large
	tooSmall.BootstrappingParameters = small
	tooSmallEval, err := NewFastEvaluator(tooSmall)
	require.NoError(t, err)
	_, _, _, err = tooSmallEval.PackAndSwitchN1ToN2([]rlwe.Ciphertext{*newFastPackingCiphertext(large, 4, 1, 523)})
	require.Error(t, err)
}
