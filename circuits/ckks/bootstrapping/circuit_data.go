package bootstrapping

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// buildBootstrapCircuitData is the pure planning and matrix-generation part
// of the Standard bootstrap constructor. Keeping it shared makes the Fast
// circuit use the same Mod1-dependent DFT scaling and level schedule without
// constructing a Standard evaluator or any evaluation-key machinery.
func buildBootstrapCircuitData(btpParams Parameters) (adjusted Parameters, mod1Params mod1.Parameters, c2s, s2c dft.Matrix, err error) {
	adjusted = btpParams
	params := btpParams.BootstrappingParameters

	if mod1Params, err = mod1.NewParametersFromLiteral(params, btpParams.Mod1ParametersLiteral); err != nil {
		return
	}

	K := mod1Params.K
	qDiff := mod1Params.QDiff
	qDiv := mod1Params.ScalingFactor().Float64() / math.Exp2(math.Round(math.Log2(float64(params.Q()[0]))))
	if qDiv > 1 {
		qDiv = 1
	}

	encoder := ckks.NewEncoder(params)
	scale := params.DefaultScale().Float64()
	offset := mod1Params.ScalingFactor().Float64() / mod1Params.MessageRatio()
	c2sScaling := new(big.Float).SetFloat64(qDiv / (K * qDiff))
	s2cScaling := new(big.Float).SetFloat64(scale / offset)

	if adjusted.CoeffsToSlotsParameters.Scaling == nil {
		adjusted.CoeffsToSlotsParameters.Scaling = c2sScaling
	} else {
		adjusted.CoeffsToSlotsParameters.Scaling = new(big.Float).Mul(adjusted.CoeffsToSlotsParameters.Scaling, c2sScaling)
	}
	if adjusted.SlotsToCoeffsParameters.Scaling == nil {
		adjusted.SlotsToCoeffsParameters.Scaling = s2cScaling
	} else {
		adjusted.SlotsToCoeffsParameters.Scaling = new(big.Float).Mul(adjusted.SlotsToCoeffsParameters.Scaling, s2cScaling)
	}

	if c2s, err = dft.NewMatrixFromLiteral(params, adjusted.CoeffsToSlotsParameters, encoder); err != nil {
		return
	}
	if s2c, err = dft.NewMatrixFromLiteral(params, adjusted.SlotsToCoeffsParameters, encoder); err != nil {
		return
	}
	return
}
