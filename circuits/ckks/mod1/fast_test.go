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

func fastMod1Q012TestParameters(t *testing.T) (ckks.Parameters, Parameters) {
	t.Helper()
	q3Generator := ring.NewNTTFriendlyPrimesGenerator(45, 2*2*16)
	q3, err := q3Generator.NextAlternatingPrimes(1)
	require.NoError(t, err)
	remainingGenerator := ring.NewNTTFriendlyPrimesGenerator(60, 2*2*16)
	remainingQ, err := remainingGenerator.NextAlternatingPrimes(6)
	require.NoError(t, err)
	q := append([]uint64{72057594037616641, 549755731969, 549756026881}, q3...)
	q = append(q, remainingQ...)
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		Q:               q,
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
	values := make([]complex128, params.MaxSlots())
	return newFastMod1EncodedInput(t, params, level, rlwe.NewScale(math.Exp2(60)), values)
}

func newFastMod1EncodedInput(t *testing.T, params ckks.Parameters, level int, scale rlwe.Scale, values []complex128) *rlwe.Ciphertext {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	pt := ckks.NewPlaintext(params, level)
	pt.IsNTT = true
	pt.Scale = scale
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

func fastMod1DefaultLikeParameters(t *testing.T) (ckks.Parameters, Parameters) {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	mod1Params, err := NewParametersFromLiteral(params, ParametersLiteral{
		LevelQ:          9,
		LogScale:        60,
		Mod1Type:        CosDiscrete,
		LogMessageRatio: 2,
		K:               16,
		Mod1Degree:      30,
		DoubleAngle:     3,
	})
	require.NoError(t, err)
	return params, mod1Params
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

func prepareFastMod1Input(t *testing.T, params ckks.Parameters, mod1Params Parameters, values []complex128) *rlwe.Ciphertext {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	plaintext := ckks.NewPlaintext(params, params.MaxLevel())
	plaintext.Scale = params.DefaultScale()
	require.NoError(t, encoder.Encode(values, plaintext))

	input := ckks.NewCiphertext(params, 1, params.MaxLevel())
	*input.MetaData = *plaintext.MetaData
	input.Value[0].Copy(plaintext.Value)
	input.Value[1].Zero()
	standard := ckks.NewEvaluator(params, nil)

	// Match the preparation used by the Standard Mod1 test: scale the
	// message to Q/MessageRatio, then to ScalingFactor/MessageRatio, apply
	// the K*QDiff normalization, and rescale to the Mod1 input level.
	scale := rlwe.NewScale(math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / mod1Params.MessageRatio()))))
	scale = scale.Div(input.Scale)
	require.NoError(t, standard.ScaleUp(input, rlwe.NewScale(math.Round(scale.Float64())), input))
	scale = mod1Params.ScalingFactor().Div(input.Scale)
	scale = scale.Div(rlwe.NewScale(mod1Params.MessageRatio()))
	require.NoError(t, standard.ScaleUp(input, rlwe.NewScale(math.Round(scale.Float64())), input))
	require.NoError(t, standard.Mul(input, 1/(mod1Params.K*mod1Params.QDiff), input))
	require.NoError(t, standard.Rescale(input, input))
	require.Equal(t, mod1Params.LevelQ, input.Level())
	return input
}

func TestFastMod1MatchesStandardCosDiscrete(t *testing.T) {
	params, mod1Params := fastMod1DefaultLikeParameters(t)
	values := make([]complex128, params.MaxSlots())
	pattern := []float64{0, 1, -2, 4, -7, 10, -13, 15}
	fraction := []float64{0.25, -0.5, 0.75, -0.125, 0.9, -0.7, 0.33, -0.2}
	q := mod1Params.QDiff * mod1Params.MessageRatio()
	for i := range values {
		values[i] = complex(pattern[i%len(pattern)]*q+fraction[i%len(fraction)], 0)
	}
	prepared := prepareFastMod1Input(t, params, mod1Params, values)
	standardInput := prepared.CopyNew()
	fastInput := prepared.CopyNew()
	for d := range fastInput.Value {
		params.RingQ().AtLevel(fastInput.Level()).MForm(fastInput.Value[d], fastInput.Value[d])
	}
	fastInput.IsMontgomery = true

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	standardEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk)))
	standard, err := NewEvaluator(standardEval, ckkspolynomial.NewEvaluator(params, standardEval), mod1Params).EvaluateNew(standardInput)
	require.NoError(t, err)

	fastEval := NewFastEvaluator(fastckks.NewEvaluator(params), ckkspolynomial.NewFastEvaluator(params, nil), mod1Params)
	fast, err := fastEval.EvaluateNew(fastInput)
	require.NoError(t, err)
	require.Equal(t, standard.Level(), fast.Level())
	require.True(t, standard.Scale.Equal(fast.Scale))
	require.True(t, fast.IsNTT)
	require.True(t, fast.IsMontgomery)
	require.Equal(t, 1, standard.Level())
	standardLow := standard.CopyNew()
	fastLow := fast.CopyNew()
	standardValues := decodeFastMod1Plaintext(t, params, standardLow)
	fastValues := decodeFastMod1Plaintext(t, params, fastLow)
	require.NotEqual(t, real(standardValues[0]), real(standardValues[1]))
	const ckksTolerance = 0.7
	for i := range fastValues {
		// A nonzero supported-domain input exercises normalization, the
		// CosDiscrete offset, polynomial evaluation, and DoubleAngle path.
		require.False(t, math.IsNaN(real(fastValues[i])))
		require.False(t, math.IsNaN(imag(fastValues[i])))
		require.InDelta(t, real(standardValues[i]), real(fastValues[i]), ckksTolerance)
		require.InDelta(t, imag(standardValues[i]), imag(fastValues[i]), ckksTolerance)
	}
}

func TestFastMod1FormalDegree30PolynomialRegression(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{56, 39, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	mod1Params, err := NewParametersFromLiteral(params, ParametersLiteral{
		LevelQ:          12,
		LogScale:        60,
		Mod1Type:        CosDiscrete,
		LogMessageRatio: 10,
		K:               16,
		Mod1Degree:      30,
		DoubleAngle:     3,
	})
	require.NoError(t, err)
	require.True(t, normalizedLogN13Profile(&params, mod1Params))

	values := make([]complex128, params.MaxSlots())
	for i := range values {
		values[i] = complex(float64((i%7)-3)*1e-5, float64((i%5)-2)*1e-5)
	}
	fastInput := newFastMod1EncodedInput(t, params, mod1Params.LevelQ, rlwe.NewScale(math.Exp2(60)), values)
	fastBefore := fastInput.CopyNew()
	standardInput := fastInput.CopyNew()
	params.RingQ().AtLevel(standardInput.Level()).IMForm(standardInput.Value[0], standardInput.Value[0])
	params.RingQ().AtLevel(standardInput.Level()).IMForm(standardInput.Value[1], standardInput.Value[1])
	standardInput.IsMontgomery = false

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	standardEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk)))
	standard, err := NewEvaluator(standardEval, ckkspolynomial.NewEvaluator(params, standardEval), mod1Params).EvaluateNew(standardInput)
	require.NoError(t, err)
	fastCKKSEval := fastckks.NewEvaluator(params)
	var capacityCheckpoints []fastckks.QPrefixCapacitySnapshot
	fastCKKSEval.SetQPrefixCapacityObserver(func(snapshot fastckks.QPrefixCapacitySnapshot) error {
		capacityCheckpoints = append(capacityCheckpoints, snapshot)
		return nil
	})
	fastPolynomialEval := ckkspolynomial.NewFastEvaluator(params, fastCKKSEval)
	fast, err := NewFastEvaluator(fastCKKSEval, fastPolynomialEval, mod1Params).EvaluateNewQPrefixRows(fastInput, fastckks.MaxQPrefixWidth)
	require.NoError(t, err)
	checkpointNames := make(map[string]bool, len(capacityCheckpoints))
	for _, checkpoint := range capacityCheckpoints {
		t.Logf("capacity name=%s L=%d rows=%d scale=%s degree=%d max_abs=%v prefix_product=%s strict_2B_lt_SQ=%t", checkpoint.Name, checkpoint.Level, checkpoint.Rows, checkpoint.Scale, checkpoint.Degree, checkpoint.MaxAbs, checkpoint.PrefixProduct, checkpoint.StrictFit)
		checkpointNames[checkpoint.Name] = true
		require.True(t, checkpoint.StrictFit, "%s exact centered component bound must satisfy strict Q-prefix capacity", checkpoint.Name)
		require.Equal(t, 4, checkpoint.Rows, "%s must retain q0123 on the accepted P93 path", checkpoint.Name)
		require.Len(t, checkpoint.MaxAbs, checkpoint.Degree+1, "%s must report every ciphertext component", checkpoint.Name)
		product, err := fastckks.QPrefixProduct(params.Q(), checkpoint.Level)
		require.NoError(t, err)
		require.Equal(t, product.String(), checkpoint.PrefixProduct)
	}
	for _, name := range []string{
		"evalmod-entry", "generated-power-2",
		"ps-baby-0", "ps-baby-1", "ps-baby-2", "ps-baby-3", "ps-baby-4",
		"ps-giant-0", "ps-giant-1", "ps-giant-2", "ps-giant-3",
		"polynomial-before-final-rescale", "polynomial-after-final-rescale",
		"double-angle-0-before", "double-angle-0-after-rescale",
		"double-angle-1-before", "double-angle-1-after-rescale",
		"double-angle-2-before", "double-angle-2-after-rescale", "evalmod-output",
	} {
		require.True(t, checkpointNames[name], "missing required P93 checkpoint %s", name)
	}
	fastCKKSEval.SetQPrefixCapacityObserver(nil)
	poisonedInput := fastInput.CopyNew()
	q3 := params.Q()[3]
	for i, residue := range poisonedInput.Value[0].Coeffs[3] {
		poisonedInput.Value[0].Coeffs[3][i] = (residue + uint64(97+i)) % q3
	}
	poisonedOutput, err := NewFastEvaluator(fastCKKSEval, fastPolynomialEval, mod1Params).EvaluateNewQPrefixRows(poisonedInput, fastckks.MaxQPrefixWidth)
	require.NoError(t, err)
	require.Equal(t, fast.Level(), poisonedOutput.Level())
	require.NotEqual(t, fast.Value[0].Coeffs[3], poisonedOutput.Value[0].Coeffs[3], "EvalMod output q3 must follow the q3 input history")
	guardEvidence := fastPolynomialEval.LastGuardSelectionEvidence()
	require.Zero(t, guardEvidence.Operations)
	require.Equal(t, -1, guardEvidence.PlanIndex)
	require.Zero(t, guardEvidence.ScalarDegree)
	require.False(t, guardEvidence.Contraction)
	require.Equal(t, standard.Level(), fast.Level())
	require.Equal(t, 4, fast.Level())
	require.Equal(t, 1, fast.Degree())
	require.True(t, standard.Scale.Equal(fast.Scale))
	require.True(t, fast.IsNTT)
	require.True(t, fast.IsMontgomery)
	outputRows, err := fastckks.QPrefixWidth(fast.Level())
	require.NoError(t, err)
	require.Equal(t, 4, outputRows)
	for component := range fast.Value {
		for row := 0; row < outputRows; row++ {
			require.Len(t, fast.Value[component].Coeffs[row], params.N())
		}
		for row := outputRows; row <= fast.Level(); row++ {
			require.Empty(t, fast.Value[component].Coeffs[row])
		}
	}
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, fastBefore.Value[d].Coeffs[limb], fastInput.Value[d].Coeffs[limb])
		}
	}

	standardValues := decodeFastMod1Plaintext(t, params, standard)
	fastValues := decodeFastMod1MaintainedPlaintext(t, params, fast)
	for i := range values {
		require.InDelta(t, real(standardValues[i]), real(fastValues[i]), 1e-2, "real slot %d", i)
		require.InDelta(t, imag(standardValues[i]), imag(fastValues[i]), 1e-2, "imag slot %d", i)
	}
}

func decodeFastMod1MaintainedPlaintext(t *testing.T, params ckks.Parameters, ct *rlwe.Ciphertext) []complex128 {
	t.Helper()
	level := 1
	source := params.RingQ().AtLevel(level).NewPoly()
	for limb := 0; limb < 2; limb++ {
		copy(source.Coeffs[limb], ct.Value[0].Coeffs[limb])
	}
	params.RingQ().AtLevel(level).INTT(source, source)
	params.RingQ().AtLevel(level).IMForm(source, source)

	q01Params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            params.LogN(),
		Q:               append([]uint64(nil), params.Q()[:2]...),
		LogDefaultScale: params.LogDefaultScale(),
	})
	require.NoError(t, err)
	encoder := ckks.NewEncoder(q01Params)
	pt := ckks.NewPlaintext(q01Params, level)
	*pt.MetaData = *ct.MetaData
	pt.Value.Copy(source)
	pt.IsNTT = false
	pt.IsMontgomery = false
	decoded := make([]complex128, pt.Slots())
	require.NoError(t, encoder.Decode(pt, decoded))
	return decoded
}

func TestFastMod1IgnoresDormantResiduesAndPreservesInput(t *testing.T) {
	params, mod1Params := fastMod1Q012TestParameters(t)
	input := newFastMod1PlaintextInput(t, params, mod1Params.LevelQ)
	require.Equal(t, 3, fastckks.MaintainedLimbCount(&params, input.Level()), "fixture must exercise legacy q012 authority")
	before := input.CopyNew()
	poisoned := input.CopyNew()
	q3 := params.Q()[3]
	for i, residue := range poisoned.Value[0].Coeffs[3] {
		poisoned.Value[0].Coeffs[3][i] = (residue + uint64(97+i)) % q3
	}

	eval := NewFastEvaluator(fastckks.NewEvaluator(params), ckkspolynomial.NewFastEvaluator(params, nil), mod1Params)
	got, err := eval.EvaluateNew(input)
	require.NoError(t, err)
	poisonedGot, err := eval.EvaluateNew(poisoned)
	require.NoError(t, err)
	for d := 0; d <= 1; d++ {
		outputRows := fastckks.MaintainedLimbCount(&params, got.Level())
		for limb := 0; limb < outputRows; limb++ {
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
	_, err = eval.EvaluateNewQPrefixRows(input, 0)
	require.Error(t, err, "an invalid explicit row count must not fall back to legacy authority")

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
