package bootstrapping

import (
	"testing"

	"github.com/stretchr/testify/require"

	ckkslintrans "github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

func fastLogN13CompressionParameters(t testing.TB) (Parameters, ckks.Parameters) {
	t.Helper()
	residual, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 39},
		LogDefaultScale: 45,
		Xs:              ring.Ternary{H: 192},
	})
	require.NoError(t, err)
	zero := 0
	logN, logSlots := 13, 12
	params, err := NewParametersFromLiteral(residual, ParametersLiteral{
		LogN:     &logN,
		LogSlots: &logSlots,
		LogP:     []int{61, 61, 61, 61, 61},
		SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{39}, {39}, {39}},
		CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{56}, {56}, {56}, {56}},
		EvalModLogScale:       pointyInt(60),
		Mod1Degree:            pointyInt(30),
		DoubleAngle:           pointyInt(3),
		K:                     pointyInt(16),
		LogMessageRatio:       pointyInt(10),
		Mod1InvDegree:         &zero,
		EphemeralSecretWeight: &zero,
	})
	require.NoError(t, err)
	params.CircuitOrder = ModUpThenEncode
	params.ResidualParameters = residual
	return params, residual
}

func pointyInt(value int) *int { return &value }

func TestFastLogN13C2SCompressionPreparation(t *testing.T) {
	params, _ := fastLogN13CompressionParameters(t)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.NoError(t, eval.ensureFastBootstrapCircuit())
	require.True(t, eval.C2SCompressionActive)
	require.Equal(t, []int{4, 2, 0, 0}, eval.C2SRestorePlan)

	originalParams := params
	_, _, original, _, err := buildBootstrapCircuitData(originalParams)
	require.NoError(t, err)
	require.Equal(t, original.MatrixLiteral, eval.C2SDFTMatrix.MatrixLiteral)
	require.Equal(t, original.Levels, eval.C2SDFTMatrix.Levels)
	require.Len(t, eval.C2SDFTMatrix.Matrices, 4)
	for i := range original.Matrices {
		got := eval.C2SDFTMatrix.Matrices[i]
		want := original.Matrices[i]
		require.Equal(t, want.LevelQ, got.LevelQ)
		require.Equal(t, want.LevelP, got.LevelP)
		require.Equal(t, want.LogDimensions, got.LogDimensions)
		require.Equal(t, want.LogBabyStepGiantStepRatio, got.LogBabyStepGiantStepRatio)
		require.Equal(t, want.N1, got.N1)
		require.ElementsMatch(t, utils.GetKeys(want.Vec), utils.GetKeys(got.Vec))
		factor := 1.0
		if i == 0 {
			factor = 16
		} else if i == 1 {
			factor = 4
		}
		require.True(t, got.Scale.Mul(rlwe.NewScale(factor)).Equal(want.Scale), "group %d scale", i)
	}
}

func TestFastLogN13C2SCompressionGuardRejectsChangedProfile(t *testing.T) {
	params, _ := fastLogN13CompressionParameters(t)
	params.CoeffsToSlotsParameters.LogSlots = 11
	_, _, matrix, _, err := buildBootstrapCircuitData(params)
	require.NoError(t, err)
	require.False(t, matchesFastLogN13C2SProfile(params, matrix))
	prepared, plan, active, err := prepareFastLogN13C2S(params, matrix)
	require.NoError(t, err)
	require.False(t, active)
	require.Nil(t, plan)
	require.Equal(t, matrix, prepared)
}

func TestFastLogN13C2SCompressionUsesSameDiagonalStructure(t *testing.T) {
	params, _ := fastLogN13CompressionParameters(t)
	_, _, original, _, err := buildBootstrapCircuitData(params)
	require.NoError(t, err)
	prepared, plan, active, err := prepareFastLogN13C2S(params, original)
	require.NoError(t, err)
	require.True(t, active)
	require.Equal(t, []int{4, 2, 0, 0}, plan)
	mathematical := original.MatrixLiteral.GenMatrices(params.BootstrappingParameters.LogN(), params.BootstrappingParameters.EncodingPrecision())
	for i := range prepared.Matrices {
		require.ElementsMatch(t, mathematical[i].DiagonalsIndexList(), utils.GetKeys(original.Matrices[i].Vec), "group %d diagonal set", i)
		require.IsType(t, ckkslintrans.LinearTransformation{}, prepared.Matrices[i])
	}
}
