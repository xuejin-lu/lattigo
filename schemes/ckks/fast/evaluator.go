package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// Evaluator is an explicit Fast CKKS execution boundary. It intentionally
// does not replace or modify ckks.Evaluator: callers select Fast execution by
// constructing this evaluator explicitly.
type Evaluator struct {
	Parameters ckks.Parameters
}

// NewEvaluator creates an explicit q0/q1-authoritative evaluator.
func NewEvaluator(params ckks.Parameters) *Evaluator {
	return &Evaluator{Parameters: params}
}

// Mul multiplies two degree-one ciphertexts and truncates the resulting
// degree-two ciphertext under Fast's zero-secret semantics. Inputs must be in
// the coefficient domain; q0 and q1 are the only authoritative limbs.
func (eval *Evaluator) Mul(op0, op1, opOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if op0 == nil || op1 == nil || opOut == nil {
		return errors.New("op0, op1 and opOut cannot be nil")
	}
	if op0.MetaData == nil || op1.MetaData == nil || opOut.MetaData == nil {
		return errors.New("op0, op1 and opOut metadata cannot be nil")
	}
	if op0.Degree() != 1 || op1.Degree() != 1 {
		return fmt.Errorf("Fast Mul requires degree-one ciphertexts")
	}
	if op0.IsNTT || op1.IsNTT {
		return errors.New("Fast Mul requires coefficient-domain ciphertexts")
	}
	if op0.IsMontgomery || op1.IsMontgomery {
		return errors.New("Fast Mul requires non-Montgomery coefficient-domain ciphertexts")
	}
	if op0.N() != eval.Parameters.N() || op1.N() != eval.Parameters.N() || opOut.N() != eval.Parameters.N() {
		return errors.New("ciphertext dimensions do not match Fast evaluator parameters")
	}
	if op0.IsBatched != op1.IsBatched {
		return errors.New("Fast Mul requires equal batching metadata")
	}

	level := utils.Min(op0.Level(), op1.Level())
	level = utils.Min(level, opOut.Level())
	if level < 1 {
		return errors.New("Fast Mul requires at least q0 and q1")
	}

	ringQ := eval.Parameters.RingQ().AtLevel(level)
	// Keep all products in level-one temporary polynomials. This allocates
	// storage only for q0/q1 and protects both inputs for aliasing cases.
	c0 := ring.NewPoly(eval.Parameters.N(), 1)
	c1a := ring.NewPoly(eval.Parameters.N(), 1)
	c1b := ring.NewPoly(eval.Parameters.N(), 1)

	if err := FastMulQ01Authoritative(ringQ, op0.Value[0], op1.Value[0], c0); err != nil {
		return fmt.Errorf("FastMulQ01Authoritative(c0): %w", err)
	}
	if err := FastMulQ01Authoritative(ringQ, op0.Value[0], op1.Value[1], c1a); err != nil {
		return fmt.Errorf("FastMulQ01Authoritative(c1a): %w", err)
	}
	if err := FastMulQ01Authoritative(ringQ, op0.Value[1], op1.Value[0], c1b); err != nil {
		return fmt.Errorf("FastMulQ01Authoritative(c1b): %w", err)
	}
	for _, limb := range []int{0, 1} {
		ringQ.SubRings[limb].Add(c1a.Coeffs[limb], c1b.Coeffs[limb], c1a.Coeffs[limb])
	}

	opOut.Resize(2, level)
	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale.Mul(op1.Scale)
	opOut.IsNTT = false
	opOut.IsMontgomery = false
	opOut.IsBatched = op0.IsBatched
	opOut.LogDimensions.Rows = utils.Max(op0.LogDimensions.Rows, op1.LogDimensions.Rows)
	opOut.LogDimensions.Cols = utils.Max(op0.LogDimensions.Cols, op1.LogDimensions.Cols)

	copyQ01(c0, opOut.Value[0])
	copyQ01(c1a, opOut.Value[1])
	return FastTruncateDegree2To1(opOut, opOut)
}

// MulNew returns a newly allocated degree-one Fast ciphertext.
func (eval *Evaluator) MulNew(op0, op1 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil || op0 == nil || op1 == nil {
		return nil, errors.New("evaluator and operands cannot be nil")
	}
	level := utils.Min(op0.Level(), op1.Level())
	opOut := ckks.NewCiphertext(eval.Parameters, 2, level)
	return opOut, eval.Mul(op0, op1, opOut)
}
