package mod1

import (
	"errors"
	"fmt"
	"math/big"

	ckkspolynomial "github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// FastEvaluator is the bounded Stage-A Fast Mod1 evaluator. It supports the
// Standard-ring CosDiscrete path without an inverse polynomial and without a
// scaling factor other than one. It deliberately does not implement
// schemes.Evaluator and never instantiates a Standard CKKS evaluator.
type FastEvaluator struct {
	Parameters          Parameters
	FastCKKS            *fastckks.Evaluator
	PolynomialEvaluator *ckkspolynomial.FastEvaluator
}

// NewFastEvaluator creates a bounded Fast Mod1 evaluator.
func NewFastEvaluator(eval *fastckks.Evaluator, evalPoly *ckkspolynomial.FastEvaluator, params Parameters) *FastEvaluator {
	return &FastEvaluator{Parameters: params, FastCKKS: eval, PolynomialEvaluator: evalPoly}
}

// EvaluateNew evaluates the supported Fast Mod1 circuit while preserving the
// Standard normalization, Chebyshev offset, target-scale schedule, and
// DoubleAngle recurrence.
func (eval *FastEvaluator) EvaluateNew(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if err := eval.validate(ct); err != nil {
		return nil, err
	}

	params := eval.FastCKKS.GetParameters()
	mod1Params := eval.Parameters
	inputScale := ct.Scale
	res := cloneFastCiphertext(*params, ct)

	// Normalize the modular reduction to mod 1 by changing the scale
	// interpretation, exactly as Standard Mod1 does.
	res.Scale = mod1Params.ScalingFactor()

	// Compute the polynomial target scale using the same Q schedule as the
	// Standard evaluator. The index is intentionally derived from the current
	// level and generated polynomial depth rather than a fixed chain layout.
	qi := params.Q()
	targetScale := res.Scale
	for i := 0; i < mod1Params.DoubleAngle; i++ {
		index := res.Level() - mod1Params.Mod1Poly.Depth() - mod1Params.DoubleAngle + i + 1
		if index < 0 || index >= len(qi) {
			return nil, fmt.Errorf("Fast Mod1 target-scale Q index %d is outside the parameter chain", index)
		}
		targetScale = targetScale.Mul(rlwe.NewScale(qi[index]))
		targetScale.Value.Sqrt(&targetScale.Value)
	}

	// Preserve the caller-side Chebyshev change of variable used by Standard
	// Mod1. The Fast scalar Add changes only the maintained q0/q1 residues.
	offset := new(big.Float).Sub(&mod1Params.Mod1Poly.B, &mod1Params.Mod1Poly.A)
	offset.Mul(offset, new(big.Float).SetFloat64(mod1Params.IntervalShrinkFactor()))
	offset.Quo(new(big.Float).SetFloat64(-0.5), offset)
	if err := eval.FastCKKS.Add(res, offset, res); err != nil {
		return nil, fmt.Errorf("Fast Mod1 cosine offset: %w", err)
	}

	polynomialResult, err := eval.PolynomialEvaluator.Evaluate(res, mod1Params.Mod1Poly, targetScale)
	if err != nil {
		return nil, fmt.Errorf("Fast Mod1 polynomial evaluation: %w", err)
	}
	// The polynomial evaluator returns an independently-owned public result;
	// keep its planned target scale through the DoubleAngle chain.
	res = polynomialResult

	sqrt2pi := mod1Params.Sqrt2Pi
	for i := 0; i < mod1Params.DoubleAngle; i++ {
		sqrt2pi *= sqrt2pi

		if err := eval.FastCKKS.MulRelin(res, res, res); err != nil {
			return nil, fmt.Errorf("Fast Mod1 double angle %d multiply: %w", i, err)
		}
		if err := eval.FastCKKS.Add(res, res, res); err != nil {
			return nil, fmt.Errorf("Fast Mod1 double angle %d doubling: %w", i, err)
		}
		if err := eval.FastCKKS.Add(res, -sqrt2pi, res); err != nil {
			return nil, fmt.Errorf("Fast Mod1 double angle %d offset: %w", i, err)
		}
		if err := eval.FastCKKS.Rescale(res, res); err != nil {
			return nil, fmt.Errorf("Fast Mod1 double angle %d rescale: %w", i, err)
		}
	}

	res.Scale = inputScale
	return res, nil
}

func (eval *FastEvaluator) validate(ct *rlwe.Ciphertext) error {
	if eval == nil || eval.FastCKKS == nil || eval.PolynomialEvaluator == nil {
		return errors.New("Fast Mod1 evaluator and dependencies cannot be nil")
	}
	if ct == nil || ct.MetaData == nil {
		return errors.New("Fast Mod1 ciphertext and metadata cannot be nil")
	}
	if eval.FastCKKS.GetParameters().RingType() != ring.Standard || eval.Parameters.Mod1Type != CosDiscrete {
		return errors.New("Fast Mod1 requires the Standard ring and Mod1Type CosDiscrete")
	}
	if eval.Parameters.Mod1InvPoly != nil {
		return errors.New("Fast Mod1 does not support Mod1 inverse polynomials")
	}
	if ct.N() != eval.FastCKKS.GetParameters().N() || ct.Degree() != 1 {
		return errors.New("Fast Mod1 requires a matching degree-one ciphertext")
	}
	if ct.Level() != eval.Parameters.LevelQ {
		return fmt.Errorf("Fast Mod1 requires input level %d, got %d", eval.Parameters.LevelQ, ct.Level())
	}
	if ct.Level() < 1 {
		return errors.New("Fast Mod1 requires input level containing q0 and q1")
	}
	if !ct.IsNTT || !ct.IsMontgomery {
		return errors.New("Fast Mod1 requires NTT-domain Montgomery input")
	}
	for d := 0; d <= 1; d++ {
		if len(ct.Value) <= d || len(ct.Value[d].Coeffs) < 2 ||
			len(ct.Value[d].Coeffs[0]) != eval.FastCKKS.GetParameters().N() ||
			len(ct.Value[d].Coeffs[1]) != eval.FastCKKS.GetParameters().N() {
			return errors.New("Fast Mod1 input has invalid q0/q1 storage")
		}
	}
	if len(eval.Parameters.Mod1Poly.Coeffs) == 0 {
		return errors.New("Fast Mod1 polynomial cannot be empty")
	}
	if eval.Parameters.Mod1Poly.Basis != bignum.Chebyshev {
		return errors.New("Fast Mod1 requires a Chebyshev polynomial")
	}
	if eval.Parameters.Mod1Poly.Depth() == 0 {
		return errors.New("Fast Mod1 requires a non-constant Mod1 polynomial")
	}
	return nil
}

// cloneFastCiphertext copies only the authoritative q0/q1 rows. The public
// ciphertext remains structurally compatible at the logical input level, but
// dormant higher rows are never read from the Fast input.
func cloneFastCiphertext(params ckks.Parameters, src *rlwe.Ciphertext) *rlwe.Ciphertext {
	dst := fastckks.NewCiphertext(params, 1, src.Level())
	*dst.MetaData = *src.MetaData
	dst.IsNTT = src.IsNTT
	dst.IsMontgomery = src.IsMontgomery
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			copy(dst.Value[d].Coeffs[limb], src.Value[d].Coeffs[limb])
		}
	}
	return dst
}
