package bootstrapping

import (
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func fastEvalModParameters(t testing.TB, logN int) (Parameters, ckks.Parameters) {
	t.Helper()
	logQ := []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60}
	if logN == 16 {
		logQ[0], logQ[1] = 54, 38
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            logQ,
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	literal := mod1.ParametersLiteral{
		LevelQ:          8,
		LogScale:        60,
		Mod1Type:        mod1.CosDiscrete,
		LogMessageRatio: 2,
		K:               4,
		Mod1Degree:      6,
		DoubleAngle:     2,
		Mod1InvDegree:   0,
	}
	return Parameters{BootstrappingParameters: params, Mod1ParametersLiteral: literal}, params
}

func fastEvalModP93Parameters(t testing.TB) (Parameters, ckks.Parameters) {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{56, 39, 40, 40, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	literal := mod1.ParametersLiteral{
		LevelQ:          12,
		LogScale:        60,
		Mod1Type:        mod1.CosDiscrete,
		LogMessageRatio: 10,
		K:               16,
		Mod1Degree:      30,
		DoubleAngle:     3,
		Mod1InvDegree:   0,
	}
	return Parameters{BootstrappingParameters: params, Mod1ParametersLiteral: literal}, params
}

func newFastEvalModInput(params ckks.Parameters, mod1Params mod1.Parameters) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, mod1Params.LevelQ)
	ct.IsNTT = true
	ct.IsMontgomery = true
	ct.Scale = mod1Params.ScalingFactor()
	return ct
}

func TestFastEvalModConstructorAndBoundary(t *testing.T) {
	bootstrapParams, params := fastEvalModParameters(t, 4)
	eval, err := NewFastEvaluator(bootstrapParams)
	require.NoError(t, err)
	require.NotNil(t, eval.Mod1Evaluator)
	require.NotNil(t, eval.PolynomialEvaluator)
	require.Equal(t, bootstrapParams.Mod1ParametersLiteral.LevelQ, eval.Mod1Parameters.LevelQ)
	require.Equal(t, bootstrapParams.Mod1ParametersLiteral.Mod1Type, eval.Mod1Parameters.Mod1Type)
	require.NotEmpty(t, eval.Mod1Parameters.Mod1Poly.Coeffs)
	require.Nil(t, eval.Mod1Parameters.Mod1InvPoly)

	input := newFastEvalModInput(params, eval.Mod1Parameters)
	out, err := eval.EvalMod(input)
	require.NoError(t, err)
	require.Equal(t, 1, out.Degree())
	require.Equal(t, eval.Mod1Parameters.LevelQ-eval.Mod1Parameters.Mod1Poly.Depth()-eval.Mod1Parameters.DoubleAngle, out.Level())
	require.True(t, out.Scale.Equal(params.DefaultScale()))
	require.True(t, out.IsNTT)
	require.True(t, out.IsMontgomery)
}

func TestFastEvalModRepeatedEvaluationDoesNotContaminateResults(t *testing.T) {
	bootstrapParams, params := fastEvalModParameters(t, 4)
	eval, err := NewFastEvaluator(bootstrapParams)
	require.NoError(t, err)
	firstInput := newFastEvalModInput(params, eval.Mod1Parameters)
	secondInput := newFastEvalModInput(params, eval.Mod1Parameters)
	secondInput.Value[0].Coeffs[0][0] = 7
	secondInput.Value[1].Coeffs[1][1] = 11
	first, err := eval.EvalMod(firstInput)
	require.NoError(t, err)
	firstCopy := first.CopyNew()
	_, err = eval.EvalMod(secondInput)
	require.NoError(t, err)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, firstCopy.Value[d].Coeffs[limb], first.Value[d].Coeffs[limb])
		}
	}
}

func TestFastEvalModDefaultLikeWrapperScale(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	bootstrapParams := Parameters{
		BootstrappingParameters: params,
		Mod1ParametersLiteral: mod1.ParametersLiteral{
			LevelQ:          9,
			LogScale:        60,
			Mod1Type:        mod1.CosDiscrete,
			LogMessageRatio: 2,
			K:               16,
			Mod1Degree:      30,
			DoubleAngle:     3,
			Mod1InvDegree:   0,
		},
	}
	eval, err := NewFastEvaluator(bootstrapParams)
	require.NoError(t, err)
	input := newFastEvalModInput(params, eval.Mod1Parameters)
	out, err := eval.EvalMod(input)
	require.NoError(t, err)
	require.Equal(t, 1, out.Degree())
	require.Equal(t, 1, out.Level())
	require.True(t, out.Scale.Equal(params.DefaultScale()))
	require.True(t, out.IsNTT)
	require.True(t, out.IsMontgomery)
}

func TestFastEvalModP93ExplicitQPrefixConsumesQ3(t *testing.T) {
	bootstrapParams, params := fastEvalModP93Parameters(t)
	eval, err := NewFastEvaluator(bootstrapParams)
	require.NoError(t, err)
	input := newFastEvalModInput(params, eval.Mod1Parameters)
	poisoned := input.CopyNew()
	q3 := params.Q()[3]
	for i, residue := range poisoned.Value[0].Coeffs[3] {
		poisoned.Value[0].Coeffs[3][i] = (residue + uint64(97+i)) % q3
	}

	cleanOutput, err := eval.EvalMod(input)
	require.NoError(t, err)
	poisonedOutput, err := eval.EvalMod(poisoned)
	require.NoError(t, err)
	require.Equal(t, 4, cleanOutput.Level())
	require.Equal(t, cleanOutput.Level(), poisonedOutput.Level())
	require.True(t, cleanOutput.Scale.Equal(params.DefaultScale()))
	q3Changed := false
	for i := range cleanOutput.Value[0].Coeffs[3] {
		if cleanOutput.Value[0].Coeffs[3][i] != poisonedOutput.Value[0].Coeffs[3][i] {
			q3Changed = true
			break
		}
	}
	require.True(t, q3Changed, "Bootstrap EvalMod explicit path must consume and propagate q3")
}
