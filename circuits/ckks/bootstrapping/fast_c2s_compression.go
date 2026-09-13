package bootstrapping

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	ckkslintrans "github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var logN13C2SRestorePlan = []int{4, 2, 0, 0}

// prepareFastLogN13C2S installs the already validated bounded capacity
// workaround. It returns the original matrix and a nil plan for every
// non-matching profile.
func prepareFastLogN13C2S(params Parameters, matrix dft.Matrix) (dft.Matrix, []int, bool, error) {
	if !matchesFastLogN13C2SProfile(params, matrix) {
		return matrix, nil, false, nil
	}

	plan := append([]int(nil), logN13C2SRestorePlan...)
	mathematical := matrix.MatrixLiteral.GenMatrices(params.BootstrappingParameters.LogN(), params.BootstrappingParameters.EncodingPrecision())
	if len(mathematical) != len(matrix.Matrices) {
		return dft.Matrix{}, nil, false, fmt.Errorf("Fast LogN13 C2S mathematical matrix count %d does not match encoded count %d", len(mathematical), len(matrix.Matrices))
	}

	encoder := ckks.NewEncoder(params.BootstrappingParameters)
	for i, k := range plan {
		if k == 0 {
			continue
		}
		original := matrix.Matrices[i]
		ltParams := ckkslintrans.Parameters{
			DiagonalsIndexList:        mathematical[i].DiagonalsIndexList(),
			LevelQ:                    original.LevelQ,
			LevelP:                    original.LevelP,
			Scale:                     original.Scale.Div(rlwe.NewScale(uint64(1) << uint(k))),
			LogDimensions:             original.LogDimensions,
			LogBabyStepGiantStepRatio: original.LogBabyStepGiantStepRatio,
		}
		replacement := ckkslintrans.NewTransformation(params.BootstrappingParameters, ltParams)
		if err := ckkslintrans.Encode(encoder, mathematical[i], replacement); err != nil {
			return dft.Matrix{}, nil, false, fmt.Errorf("cannot encode compressed Fast LogN13 C2S group %d: %w", i, err)
		}
		matrix.Matrices[i] = replacement
	}
	return matrix, plan, true, nil
}

func matchesFastLogN13C2SProfile(params Parameters, matrix dft.Matrix) bool {
	bootstrap := params.BootstrappingParameters
	if bootstrap.LogN() != 13 || bootstrap.LogMaxSlots() != 12 || bootstrap.LevelsConsumedPerRescaling() != 1 {
		return false
	}
	if params.CoeffsToSlotsParameters.Type != dft.HomomorphicEncode ||
		params.CoeffsToSlotsParameters.Format != dft.RepackImagAsReal ||
		params.CoeffsToSlotsParameters.LogSlots != 12 ||
		params.CoeffsToSlotsParameters.Levels == nil || len(params.CoeffsToSlotsParameters.Levels) != 4 {
		return false
	}
	for _, factors := range params.CoeffsToSlotsParameters.Levels {
		if factors != 1 {
			return false
		}
	}
	if matrix.Type != dft.HomomorphicEncode || matrix.Format != dft.RepackImagAsReal ||
		matrix.LogSlots != 12 || len(matrix.Levels) != 4 || len(matrix.Matrices) != 4 {
		return false
	}
	for _, factors := range matrix.Levels {
		if factors != 1 {
			return false
		}
	}
	if matrix.LevelQ != 16 || params.Mod1ParametersLiteral.LevelQ != 12 || params.Mod1ParametersLiteral.Depth() != 8 {
		return false
	}
	q := bootstrap.Q()
	if len(q) != 17 {
		return false
	}
	expectedBits := []int{55, 39, 39, 39, 39, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56}
	for i, qi := range q {
		bitLen := bits.Len64(qi)
		if bitLen < expectedBits[i] || bitLen > expectedBits[i]+1 {
			return false
		}
	}
	for _, matrixPart := range matrix.Matrices {
		if matrixPart.LevelQ != matrix.LevelQ || matrixPart.LevelP != matrix.LevelP || matrixPart.LogDimensions.Cols != 12 {
			return false
		}
	}
	return true
}
