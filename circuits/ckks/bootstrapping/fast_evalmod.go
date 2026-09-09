package bootstrapping

import (
	"errors"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// CoeffsToSlots applies the Fast DFT through the ordinary bootstrap stage
// boundary. The compatibility constructor initializes the matrices lazily,
// just as the complete Fast Bootstrap path does.
func (eval *FastEvaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	if err = eval.ensureFastBootstrapCircuit(); err != nil {
		return nil, nil, err
	}
	return eval.DFTEvaluator.CoeffsToSlotsNew(ctIn, eval.C2SDFTMatrix)
}

// SlotsToCoeffs applies the Fast DFT through the ordinary bootstrap stage
// boundary.
func (eval *FastEvaluator) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if err := eval.ensureFastBootstrapCircuit(); err != nil {
		return nil, err
	}
	return eval.DFTEvaluator.SlotsToCoeffsNew(ctReal, ctImag, eval.S2CDFTMatrix)
}

// EvalMod applies the bounded Fast Mod1 circuit and restores the Bootstrap
// boundary's public default scale metadata, as Standard Bootstrap does.
func (eval *FastEvaluator) EvalMod(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil || eval.Mod1Evaluator == nil {
		return nil, errors.New("Fast Bootstrap Mod1 evaluator cannot be nil")
	}
	ctOut, err := eval.Mod1Evaluator.EvaluateNew(ctIn)
	if err != nil {
		return nil, err
	}
	ctOut.Scale = eval.Parameters.BootstrappingParameters.DefaultScale()
	return ctOut, nil
}
