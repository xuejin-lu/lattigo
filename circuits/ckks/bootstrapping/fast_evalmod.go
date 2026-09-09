package bootstrapping

import (
	"errors"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

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
