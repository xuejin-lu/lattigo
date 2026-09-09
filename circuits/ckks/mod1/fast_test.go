package mod1

import (
	"math"
	"testing"

	"github.com/stretchr/testify/require"

	ckkspolynomial "github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

func fastMod1TestParameters(t *testing.T) (ckks.Parameters, Parameters) {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	mod1Params, err := NewParametersFromLiteral(params, ParametersLiteral{
		LevelQ:          8,
		LogScale:        60,
		Mod1Type:        CosDiscrete,
		LogMessageRatio: 2,
		K:               4,
		Mod1Degree:      6,
		DoubleAngle:     2,
	})
	require.NoError(t, err)
	return params, mod1Params
}

func newFastMod1PlaintextInput(t *testing.T, params ckks.Parameters, level int) *rlwe.Ciphertext {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	pt := ckks.NewPlaintext(params, level)
	pt.IsNTT = true
	pt.Scale = rlwe.NewScale(math.Exp2(60))
	values := make([]complex128, params.MaxSlots())
	for i := range values {
		values[i] = 0
	}
	require.NoError(t, encoder.Encode(values, pt))
	for limb := 0; limb <= level; limb++ {
		params.RingQ().SubRings[limb].MForm(pt.Value.Coeffs[limb], pt.Value.Coeffs[limb])
	}
	pt.IsMontgomery = true

	ct := ckks.NewCiphertext(params, 1, level)
	*ct.MetaData = *pt.MetaData
	ct.Value[0].Copy(pt.Value)
	ct.Value[1].Zero()
	return ct
}

func decodeFastMod1Plaintext(t *testing.T, params ckks.Parameters, ct *rlwe.Ciphertext) []complex128 {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	pt := ckks.NewPlaintext(params, ct.Level())
	*pt.MetaData = *ct.MetaData
	pt.Value.Copy(ct.Value[0])
	if pt.IsMontgomery {
		params.RingQ().AtLevel(pt.Level()).INTT(pt.Value, pt.Value)
		params.RingQ().AtLevel(pt.Level()).IMForm(pt.Value, pt.Value)
		pt.IsNTT = false
		pt.IsMontgomery = false
	}
	values := make([]complex128, pt.Slots())
	require.NoError(t, encoder.Decode(pt, values))
	return values
}

func TestFastMod1MatchesStandardCosDiscrete(t *testing.T) {
	params, mod1Params := fastMod1TestParameters(t)
	input := newFastMod1PlaintextInput(t, params, mod1Params.LevelQ)

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	standardEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk)))
	standard, err := NewEvaluator(standardEval, ckkspolynomial.NewEvaluator(params, standardEval), mod1Params).EvaluateNew(input.CopyNew())
	require.NoError(t, err)

	fastEval := NewFastEvaluator(fastckks.NewEvaluator(params), ckkspolynomial.NewFastEvaluator(params, nil), mod1Params)
	fast, err := fastEval.EvaluateNew(input)
	require.NoError(t, err)
	require.Equal(t, standard.Level(), fast.Level())
	require.True(t, standard.Scale.Equal(fast.Scale))
	require.True(t, fast.IsNTT)
	require.True(t, fast.IsMontgomery)
	standardLow := standard.CopyNew()
	standardLow.Resize(standardLow.Degree(), 1)
	fastLow := fast.CopyNew()
	fastLow.Resize(fastLow.Degree(), 1)
	fastValues := decodeFastMod1Plaintext(t, params, fastLow)
	for i := range fastValues {
		// The zero-message clear oracle exercises the complete supported
		// normalization, CosDiscrete offset, polynomial, and DoubleAngle path.
		require.InDelta(t, 0, real(fastValues[i]), 1e-5)
		require.InDelta(t, 0, imag(fastValues[i]), 1e-5)
	}
	_ = standardLow
}

func TestFastMod1IgnoresDormantResiduesAndPreservesInput(t *testing.T) {
	params, mod1Params := fastMod1TestParameters(t)
	input := newFastMod1PlaintextInput(t, params, mod1Params.LevelQ)
	before := input.CopyNew()
	poisoned := input.CopyNew()
	for d := range poisoned.Value {
		for limb := 2; limb <= poisoned.Level(); limb++ {
			for i := range poisoned.Value[d].Coeffs[limb] {
				poisoned.Value[d].Coeffs[limb][i] = ^uint64(0) - uint64(i+limb+d)
			}
		}
	}

	eval := NewFastEvaluator(fastckks.NewEvaluator(params), ckkspolynomial.NewFastEvaluator(params, nil), mod1Params)
	got, err := eval.EvaluateNew(input)
	require.NoError(t, err)
	poisonedGot, err := eval.EvaluateNew(poisoned)
	require.NoError(t, err)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, got.Value[d].Coeffs[limb], poisonedGot.Value[d].Coeffs[limb])
		}
	}
	require.Equal(t, got.Level(), poisonedGot.Level())
	require.True(t, got.Scale.Equal(poisonedGot.Scale))
	for d := range input.Value {
		for limb := 0; limb <= input.Level(); limb++ {
			require.Equal(t, before.Value[d].Coeffs[limb], input.Value[d].Coeffs[limb])
		}
	}
}

func TestFastMod1Validation(t *testing.T) {
	params, mod1Params := fastMod1TestParameters(t)
	eval := NewFastEvaluator(fastckks.NewEvaluator(params), ckkspolynomial.NewFastEvaluator(params, nil), mod1Params)
	input := newFastMod1PlaintextInput(t, params, mod1Params.LevelQ)

	_, err := eval.EvaluateNew(nil)
	require.Error(t, err)

	nonNTT := input.CopyNew()
	nonNTT.IsNTT = false
	_, err = eval.EvaluateNew(nonNTT)
	require.Error(t, err)

	nonMontgomery := input.CopyNew()
	nonMontgomery.IsMontgomery = false
	_, err = eval.EvaluateNew(nonMontgomery)
	require.Error(t, err)

	wrongLevel := input.CopyNew()
	wrongLevel.Resize(1, mod1Params.LevelQ-1)
	_, err = eval.EvaluateNew(wrongLevel)
	require.Error(t, err)

	degreeTwo := input.CopyNew()
	degreeTwo.Resize(2, degreeTwo.Level())
	_, err = eval.EvaluateNew(degreeTwo)
	require.Error(t, err)

	wrongRing, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 39, 39, 39, 39, 39, 39, 39, 39},
		LogDefaultScale: 30,
		RingType:        ring.ConjugateInvariant,
	})
	require.NoError(t, err)
	wrongRingEval := NewFastEvaluator(fastckks.NewEvaluator(wrongRing), ckkspolynomial.NewFastEvaluator(wrongRing, nil), mod1Params)
	wrongRingInput := newFastMod1PlaintextInput(t, wrongRing, mod1Params.LevelQ)
	_, err = wrongRingEval.EvaluateNew(wrongRingInput)
	require.Error(t, err)

	constant := mod1Params
	constant.Mod1Poly = mod1Params.Mod1Poly.Clone()
	constant.Mod1Poly.Coeffs = constant.Mod1Poly.Coeffs[:1]
	_, err = NewFastEvaluator(fastckks.NewEvaluator(params), ckkspolynomial.NewFastEvaluator(params, nil), constant).EvaluateNew(input)
	require.Error(t, err)
}
