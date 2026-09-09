package bootstrapping

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	ckkspolynomial "github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

// FastEvaluator is the explicit Stage-A Fast Bootstrap boundary. It does not
// embed a Standard CKKS evaluator or accept Bootstrap evaluation keys.
type FastEvaluator struct {
	Parameters          Parameters
	FastCKKS            *fastckks.Evaluator
	Mod1Parameters      mod1.Parameters
	PolynomialEvaluator *ckkspolynomial.FastEvaluator
	Mod1Evaluator       *mod1.FastEvaluator
	DFTEvaluator        *dft.FastEvaluator
	C2SDFTMatrix        dft.Matrix
	S2CDFTMatrix        dft.Matrix

	fastPackingInitialized   bool
	fastPackingErr           error
	fastBootstrapInitialized bool
	fastBootstrapErr         error
	xPow2N1                  []ring.Poly
	xPow2N2                  []ring.Poly
	xPow2InvN1               []ring.Poly
	xPow2InvN2               []ring.Poly
}

// NewFastEvaluator creates the Fast ScaleDown boundary without constructing
// Standard circuit evaluators or requiring evaluation keys.
func NewFastEvaluator(params Parameters) (*FastEvaluator, error) {
	if params.BootstrappingParameters.RingType() != ring.Standard {
		return nil, fmt.Errorf("Fast Bootstrap ScaleDown requires the Standard ring")
	}
	fastEval := fastckks.NewEvaluator(params.BootstrappingParameters)
	literal := params.Mod1ParametersLiteral
	var mod1Params mod1.Parameters
	if literal.K != 0 || literal.Mod1Degree != 0 || literal.DoubleAngle != 0 || literal.Mod1InvDegree != 0 {
		var err error
		mod1Params, err = mod1.NewParametersFromLiteral(params.BootstrappingParameters, literal)
		if err != nil {
			return nil, fmt.Errorf("cannot construct Fast Mod1 parameters: %w", err)
		}
	} else {
		// The existing ScaleDown-only boundary is also used with a deliberately
		// minimal literal that has no EvalMod polynomial configuration. Preserve
		// that narrow construction for callers that do not request EvalMod.
		mod1Params = mod1.Parameters{LogDefaultScale: literal.LogScale, LogMessageRatio: literal.LogMessageRatio}
	}
	polyEval := ckkspolynomial.NewFastEvaluator(params.BootstrappingParameters, fastEval)
	var mod1Eval *mod1.FastEvaluator
	if len(mod1Params.Mod1Poly.Coeffs) != 0 {
		mod1Eval = mod1.NewFastEvaluator(fastEval, polyEval, mod1Params)
	}
	return &FastEvaluator{
		Parameters:          params,
		FastCKKS:            fastEval,
		Mod1Parameters:      mod1Params,
		PolynomialEvaluator: polyEval,
		Mod1Evaluator:       mod1Eval,
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
