package bootstrapping

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

// FastEvaluator is the narrow Fast Bootstrap boundary. It deliberately does
// not embed a Standard CKKS evaluator or accept Bootstrap evaluation keys.
// ScaleDown and the Stage-A ModUp boundary are implemented here; full
// Bootstrap orchestration remains a future stage.
type FastEvaluator struct {
	Parameters     Parameters
	FastCKKS       *fastckks.Evaluator
	Mod1Parameters mod1.Parameters
}

// NewFastEvaluator creates the Fast ScaleDown boundary without constructing
// Standard circuit evaluators or requiring evaluation keys.
func NewFastEvaluator(params Parameters) (*FastEvaluator, error) {
	if params.BootstrappingParameters.RingType() != ring.Standard {
		return nil, fmt.Errorf("Fast Bootstrap ScaleDown requires the Standard ring")
	}
	return &FastEvaluator{
		Parameters: params,
		FastCKKS:   fastckks.NewEvaluator(params.BootstrappingParameters),
		Mod1Parameters: mod1.Parameters{
			LogDefaultScale: params.Mod1ParametersLiteral.LogScale,
			LogMessageRatio: params.Mod1ParametersLiteral.LogMessageRatio,
		},
	}, nil
}

// ScaleDown follows Standard Bootstrap ScaleDown's cheap DropLevel, integer
// scale alignment, and final RescaleTo sequence. It always returns Level 0.
func (eval *FastEvaluator) ScaleDown(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, *rlwe.Scale, error) {
	if err := eval.validateScaleDown(ctIn); err != nil {
		return nil, nil, err
	}
	params := &eval.Parameters.BootstrappingParameters
	r := params.RingQ()

	for ctIn.Level() != 0 && checkMessageRatio(ctIn, eval.Mod1Parameters.MessageRatio(), r) {
		ctIn.Resize(ctIn.Degree(), ctIn.Level()-1)
	}

	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[ctIn.Level()])
	currentMessageRatio = currentMessageRatio.Div(ctIn.Scale)
	targetMessageRatio := rlwe.NewScale(eval.Mod1Parameters.MessageRatio())
	scaleUp := currentMessageRatio.Div(targetMessageRatio)
	if scaleUp.Cmp(rlwe.NewScale(0.5)) == -1 {
		return nil, nil, fmt.Errorf("initial Q/Scale = %f < 0.5*Q[0]/MessageRatio = %f", currentMessageRatio.Float64(), targetMessageRatio.Float64())
	}

	scaleUpBigint := scaleUp.BigInt()
	if err := eval.FastCKKS.MulIntegerMaintained(ctIn, scaleUpBigint, ctIn); err != nil {
		return nil, nil, err
	}
	ctIn.Scale = ctIn.Scale.Mul(rlwe.NewScale(scaleUpBigint))

	targetScale := new(big.Float).SetPrec(256).SetInt(r.ModulusAtLevel[0])
	targetScale.Quo(targetScale, new(big.Float).SetFloat64(eval.Mod1Parameters.MessageRatio()))
	if ctIn.Level() != 0 {
		if err := eval.FastCKKS.RescaleTo(ctIn, rlwe.NewScale(targetScale), ctIn); err != nil {
			return nil, nil, err
		}
	}
	errScale := ctIn.Scale.Div(rlwe.NewScale(targetScale))
	return ctIn, &errScale, nil
}

func (eval *FastEvaluator) validateScaleDown(ct *rlwe.Ciphertext) error {
	if eval == nil || eval.FastCKKS == nil {
		return errors.New("Fast Bootstrap evaluator cannot be nil")
	}
	if ct == nil || ct.MetaData == nil {
		return errors.New("Fast Bootstrap ScaleDown ciphertext and metadata cannot be nil")
	}
	if eval.Parameters.BootstrappingParameters.RingType() != ring.Standard {
		return errors.New("Fast Bootstrap ScaleDown requires the Standard ring")
	}
	if ct.N() != eval.Parameters.BootstrappingParameters.N() || ct.Degree() != 1 {
		return errors.New("Fast Bootstrap ScaleDown requires a matching degree-one ciphertext")
	}
	if ct.Level() < 0 || ct.Level() > eval.Parameters.BootstrappingParameters.MaxLevel() {
		return errors.New("Fast Bootstrap ScaleDown ciphertext level is outside Bootstrap parameters")
	}
	if !ct.IsNTT {
		return errors.New("Fast Bootstrap ScaleDown requires an NTT-domain ciphertext")
	}
	if len(ct.Value) != 2 || len(ct.Value[0].Coeffs) <= ct.Level() || len(ct.Value[1].Coeffs) <= ct.Level() {
		return errors.New("Fast Bootstrap ScaleDown ciphertext storage is invalid")
	}
	return nil
}
