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
	if eval == nil || eval.FastCKKS == nil || eval.PolynomialEvaluator == nil {
		return nil, errors.New("Fast Mod1 evaluator and dependencies cannot be nil")
	}
	if ct == nil {
		return nil, errors.New("Fast Mod1 ciphertext and metadata cannot be nil")
	}
	rows := fastckks.MaintainedLimbCount(eval.FastCKKS.GetParameters(), ct.Level())
	return eval.evaluateNew(ct, rows)
}

// EvaluateNewQPrefixRows evaluates the Fast Mod1 circuit with an explicit
// arithmetic Q-prefix row count. EvaluateNew retains the legacy maintained-row
// contract for compatibility.
func (eval *FastEvaluator) EvaluateNewQPrefixRows(ct *rlwe.Ciphertext, rows int) (*rlwe.Ciphertext, error) {
	return eval.evaluateNew(ct, rows)
}

func (eval *FastEvaluator) evaluateNew(ct *rlwe.Ciphertext, authorityRows int) (*rlwe.Ciphertext, error) {
	if err := eval.validate(ct, authorityRows); err != nil {
		return nil, err
	}

	params := eval.FastCKKS.GetParameters()
	mod1Params := eval.Parameters
	inputScale := ct.Scale
	res := cloneFastCiphertext(*params, ct, authorityRows)
	rows, err := qPrefixRowsAtLevel(res.Level(), authorityRows)
	if err != nil {
		return nil, fmt.Errorf("Fast Mod1 input Q-prefix: %w", err)
	}
	if err := eval.FastCKKS.ObserveQPrefixCapacity("evalmod-entry", res, rows); err != nil {
		return nil, err
	}

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
	if err := eval.FastCKKS.AddScalarQPrefixRows(res, offset, rows, res); err != nil {
		return nil, fmt.Errorf("Fast Mod1 cosine offset: %w", err)
	}

	polynomialResult, err := eval.PolynomialEvaluator.EvaluateQPrefixRows(res, mod1Params.Mod1Poly, targetScale, authorityRows)
	if err != nil {
		return nil, fmt.Errorf("Fast Mod1 polynomial evaluation: %w", err)
	}
	// The polynomial evaluator returns an independently-owned public result;
	// keep its planned target scale through the DoubleAngle chain.
	res = polynomialResult
	rows, err = qPrefixRowsAtLevel(res.Level(), authorityRows)
	if err != nil {
		return nil, err
	}
	if err := eval.FastCKKS.ObserveQPrefixCapacity("evalmod-polynomial-output", res, rows); err != nil {
		return nil, err
	}

	sqrt2pi := mod1Params.Sqrt2Pi
	for i := 0; i < mod1Params.DoubleAngle; i++ {
		sqrt2pi *= sqrt2pi
		rows, err = qPrefixRowsAtLevel(res.Level(), authorityRows)
		if err != nil {
			return nil, err
		}
		if err := eval.FastCKKS.ObserveQPrefixCapacity(fmt.Sprintf("double-angle-%d-before", i), res, rows); err != nil {
			return nil, err
		}

		if err := eval.FastCKKS.MulRelinElementQPrefixRows(res, res.El(), rows, res); err != nil {
			return nil, fmt.Errorf("Fast Mod1 double angle %d multiply: %w", i, err)
		}
		if err := eval.FastCKKS.AddQPrefixRows(res, res, res, rows); err != nil {
			return nil, fmt.Errorf("Fast Mod1 double angle %d doubling: %w", i, err)
		}
		if err := eval.FastCKKS.AddScalarQPrefixRows(res, -sqrt2pi, rows, res); err != nil {
			return nil, fmt.Errorf("Fast Mod1 double angle %d offset: %w", i, err)
		}
		if err := eval.FastCKKS.RescaleQPrefixRows(res, rows, res); err != nil {
			return nil, fmt.Errorf("Fast Mod1 double angle %d rescale: %w", i, err)
		}
		rows, err = qPrefixRowsAtLevel(res.Level(), authorityRows)
		if err != nil {
			return nil, err
		}
		if err := eval.FastCKKS.ObserveQPrefixCapacity(fmt.Sprintf("double-angle-%d-after-rescale", i), res, rows); err != nil {
			return nil, err
		}
	}

	res.Scale = inputScale
	if err := eval.FastCKKS.ObserveQPrefixCapacity("evalmod-output", res, rows); err != nil {
		return nil, err
	}
	return res, nil
}

func qPrefixRowsAtLevel(level, authorityRows int) (int, error) {
	width, err := fastckks.QPrefixWidth(level)
	if err != nil {
		return 0, err
	}
	if authorityRows < 1 {
		return 0, fmt.Errorf("Fast Mod1 Q-prefix rows must be positive, got %d", authorityRows)
	}
	if authorityRows > width {
		authorityRows = width
	}
	return authorityRows, nil
}

func (eval *FastEvaluator) validate(ct *rlwe.Ciphertext, authorityRows int) error {
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
	width, err := fastckks.QPrefixWidth(ct.Level())
	if err != nil {
		return fmt.Errorf("Fast Mod1 input Q-prefix: %w", err)
	}
	if authorityRows < 1 || authorityRows > width {
		return fmt.Errorf("Fast Mod1 input Q-prefix rows %d exceed Level %d width %d", authorityRows, ct.Level(), width)
	}
	for d := 0; d <= 1; d++ {
		if len(ct.Value) <= d || len(ct.Value[d].Coeffs) < authorityRows {
			return errors.New("Fast Mod1 input has invalid maintained storage")
		}
		for row := 0; row < authorityRows; row++ {
			if len(ct.Value[d].Coeffs[row]) != eval.FastCKKS.GetParameters().N() {
				return fmt.Errorf("Fast Mod1 input row %d has invalid maintained storage", row)
			}
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

// cloneFastCiphertext copies only the authoritative Q-prefix rows. The public
// ciphertext remains structurally compatible at the logical input level, but
// dormant higher rows are never read from the Fast input.
func cloneFastCiphertext(params ckks.Parameters, src *rlwe.Ciphertext, rows int) *rlwe.Ciphertext {
	dst := fastckks.NewCiphertext(params, 1, src.Level())
	*dst.MetaData = *src.MetaData
	dst.IsNTT = src.IsNTT
	dst.IsMontgomery = src.IsMontgomery
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < rows; limb++ {
			copy(dst.Value[d].Coeffs[limb], src.Value[d].Coeffs[limb])
		}
	}
	return dst
}
