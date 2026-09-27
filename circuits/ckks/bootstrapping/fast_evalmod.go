package bootstrapping

import (
	"errors"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

// CoeffsToSlots applies the Fast DFT through the ordinary bootstrap stage
// boundary. The compatibility constructor initializes the matrices lazily,
// just as the complete Fast Bootstrap path does.
func (eval *FastEvaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	if err = eval.ensureFastBootstrapCircuit(); err != nil {
		return nil, nil, err
	}
	return eval.DFTEvaluator.CoeffsToSlotsNewWithRestorePlan(ctIn, eval.C2SDFTMatrix, eval.C2SRestorePlan)
}

// SlotsToCoeffs applies the Fast DFT through the ordinary bootstrap stage
// boundary.
func (eval *FastEvaluator) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil || ctReal == nil {
		return nil, errors.New("Fast SlotsToCoeffs evaluator and real ciphertext cannot be nil")
	}
	if err := eval.ensureFastBootstrapCircuit(); err != nil {
		return nil, err
	}
	rows := fastckks.MaintainedLimbCount(&eval.Parameters.BootstrappingParameters, ctReal.Level())
	return eval.DFTEvaluator.SlotsToCoeffsNewQPrefixRows(ctReal, ctImag, eval.S2CDFTMatrix, rows)
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
