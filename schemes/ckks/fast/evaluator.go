package fast

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

// Evaluator is an explicit Fast CKKS execution boundary. It intentionally
// does not replace or modify ckks.Evaluator: callers select Fast execution by
// constructing this evaluator explicitly. Evaluators reuse internal scratch
// and are intended for one execution stream at a time; use separate instances
// for concurrent evaluation.
type Evaluator struct {
	Parameters              ckks.Parameters
	rescaleScratch          fastRescaleScratch
	nttScratch              [3]ring.Poly
	linearTransformScratch  fastLinearTransformScratch
	automorphismCore        fastcore.AutomorphismCore
	addSubCore              fastcore.AddSubCore
	lastBSGSBabyRotations   int
	qPrefixCapacityObserver func(QPrefixCapacitySnapshot) error
}

// NewEvaluator creates an explicit Fast Q-prefix evaluator. Individual
// production operations continue to select the row authority defined by
// their current producer contract.
func NewEvaluator(params ckks.Parameters) *Evaluator {
	prefixWidth := qPrefixWidthOrPanic(params.MaxLevel())
	eval := &Evaluator{
		Parameters:             params,
		rescaleScratch:         newFastRescaleScratch(params.RingQ()),
		linearTransformScratch: newFastLinearTransformScratch(params.N(), prefixWidth),
		automorphismCore:       fastcore.NewAutomorphismWorkspace(params.N(), prefixWidth),
		addSubCore:             fastcore.NewAddSubWorkspace(),
	}
	for i := range eval.nttScratch {
		eval.nttScratch[i] = ring.NewPoly(params.N(), prefixWidth-1)
	}
	return eval
}

// GetParameters returns the Fast evaluator's CKKS parameters.
func (eval *Evaluator) GetParameters() *ckks.Parameters { return &eval.Parameters }

// GetRLWEParameters satisfies rlwe.ParameterProvider without exposing a
// Standard evaluator or any full-Q/QP key-switching operation.
func (eval *Evaluator) GetRLWEParameters() *rlwe.Parameters { return &eval.Parameters.Parameters }

// MulIntegerMaintained multiplies the maintained residues by an integer
// without changing the ciphertext scale. It is intentionally separate from
// Mul: Bootstrap ScaleDown must support the valid Level-0/r0-only state while
// the generic Fast arithmetic surface continues to require q0 and q1.
func (eval *Evaluator) MulIntegerMaintained(op0 *rlwe.Ciphertext, scalar *big.Int, opOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast integer multiplication evaluator, operands and scalar cannot be nil")
	}
	level := 0
	if op0 != nil {
		level = op0.Level()
	}
	if opOut != nil && opOut.Level() < level {
		level = opOut.Level()
	}
	rows := maintainedLimbCount(&eval.Parameters, level)
	return eval.mulIntegerRows(op0, scalar, opOut, rows)
}

// MulIntegerQPrefixRows multiplies exactly rows explicit Q-prefix residues by
// scalar without changing the ciphertext scale. Rows above the requested
// prefix are neither read nor written. It is intended for boundaries, such as
// Level-0 ModUp, whose row authority is wider than the legacy arithmetic path.
func (eval *Evaluator) MulIntegerQPrefixRows(op0 *rlwe.Ciphertext, scalar *big.Int, rows int, opOut *rlwe.Ciphertext) error {
	if eval == nil || op0 == nil || opOut == nil {
		return errors.New("Fast integer multiplication evaluator and operands cannot be nil")
	}
	level := min(op0.Level(), opOut.Level())
	if err := eval.validateExplicitRows(level, rows); err != nil {
		return err
	}
	if opOut.Level() == op0.Level() && opOut.Degree() == op0.Degree() {
		return eval.mulIntegerRowsPreservingHigherRows(op0, scalar, opOut, rows)
	}
	return eval.mulIntegerRows(op0, scalar, opOut, rows)
}

// MulQPrefixRows multiplies a ciphertext by a scalar over exactly rows
// authoritative Q-prefix limbs. This surface is intentionally scalar-only;
// ciphertext products remain on the separately specified multiplication
// path.
func (eval *Evaluator) MulQPrefixRows(op0 *rlwe.Ciphertext, scalar rlwe.Operand, rows int, opOut *rlwe.Ciphertext) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if op0.Level() != opOut.Level() {
		return errors.New("Fast scalar Q-prefix multiplication requires matching input/output levels")
	}
	width, err := QPrefixWidth(op0.Level())
	if err != nil {
		return err
	}
	if rows < 1 || rows > width {
		return fmt.Errorf("Fast scalar multiplication row count %d must be in [1,%d] at level %d", rows, width, op0.Level())
	}
	c, err := eval.scalar(scalar)
	if err != nil {
		return err
	}
	scale := rlwe.NewScale(1)
	if !c.IsInt() {
		scale, err = eval.coefficientScale(op0.Level())
		if err != nil {
			return err
		}
	}
	return eval.mulScalarAtScaleRows(op0, c, scale, opOut, rows)
}

// MulElementQPrefixRows multiplies an RLWE element operand using exactly rows
// authoritative Q-prefix residues.
func (eval *Evaluator) MulElementQPrefixRows(op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], rows int, opOut *rlwe.Ciphertext) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if err := eval.validateBinaryPrefix(op0, op1, opOut); err != nil {
		return err
	}
	level := min(op0.Level(), min(op1.Level(), opOut.Level()))
	if err := eval.validateExplicitRows(level, rows); err != nil {
		return err
	}
	return eval.mulElementRows(op0, op1.El(), opOut, false, rows)
}

// MulRelinElementQPrefixRows multiplies and applies Fast's zero-secret
// degree-two truncation over exactly rows authoritative Q-prefix residues.
func (eval *Evaluator) MulRelinElementQPrefixRows(op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], rows int, opOut *rlwe.Ciphertext) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if err := eval.validateBinaryPrefix(op0, op1, opOut); err != nil {
		return err
	}
	level := min(op0.Level(), min(op1.Level(), opOut.Level()))
	if err := eval.validateExplicitRows(level, rows); err != nil {
		return err
	}
	return eval.mulElementRows(op0, op1.El(), opOut, true, rows)
}

// RelinearizeQPrefixRows truncates c2 while copying exactly rows residues.
func (eval *Evaluator) RelinearizeQPrefixRows(op0, opOut *rlwe.Ciphertext, rows int) error {
	if err := eval.validateUnary(op0, opOut); err != nil {
		return err
	}
	if err := eval.validateExplicitRows(op0.Level(), rows); err != nil {
		return err
	}
	return fastTruncateDegree2To1Rows(eval.Parameters.RingQ(), op0, opOut, rows)
}

func (eval *Evaluator) validateExplicitRows(level, rows int) error {
	width, err := QPrefixWidth(level)
	if err != nil {
		return err
	}
	if rows < 1 || rows > width || rows > level+1 {
		return fmt.Errorf("Fast explicit row count %d must be in [1,%d] at level %d", rows, min(width, level+1), level)
	}
	return nil
}

func (eval *Evaluator) mulIntegerRows(op0 *rlwe.Ciphertext, scalar *big.Int, opOut *rlwe.Ciphertext, rows int) error {
	return eval.mulIntegerRowsWithPolicy(op0, scalar, opOut, rows, false)
}

func (eval *Evaluator) mulIntegerRowsPreservingHigherRows(op0 *rlwe.Ciphertext, scalar *big.Int, opOut *rlwe.Ciphertext, rows int) error {
	return eval.mulIntegerRowsWithPolicy(op0, scalar, opOut, rows, true)
}

func (eval *Evaluator) mulIntegerRowsWithPolicy(op0 *rlwe.Ciphertext, scalar *big.Int, opOut *rlwe.Ciphertext, rows int, preserveHigherRows bool) error {
	if eval == nil || op0 == nil || opOut == nil || scalar == nil {
		return errors.New("Fast integer multiplication evaluator, operands and scalar cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return errors.New("Fast integer multiplication metadata cannot be nil")
	}
	if op0.N() != eval.Parameters.N() || opOut.N() != eval.Parameters.N() {
		return errors.New("Fast integer multiplication dimensions do not match parameters")
	}
	if !op0.IsNTT || !opOut.IsNTT {
		return errors.New("Fast integer multiplication requires NTT-domain operands")
	}
	if op0.IsMontgomery != opOut.IsMontgomery {
		return errors.New("Fast integer multiplication requires matching Montgomery representations")
	}

	level := op0.Level()
	if opOut.Level() < level {
		level = opOut.Level()
	}
	ringQ := eval.Parameters.RingQ()
	for d := range op0.Value {
		if err := validatePrefixRows(ringQ, level, rows, op0.Value[d]); err != nil {
			return fmt.Errorf("input component %d: %w", d, err)
		}
	}
	if !preserveHigherRows {
		Resize(opOut, op0.Degree(), level, eval.Parameters.N())
	}
	for d := range opOut.Value {
		if err := validatePrefixRows(ringQ, level, rows, opOut.Value[d]); err != nil {
			return fmt.Errorf("output component %d: %w", d, err)
		}
	}
	var scalarMontgomery [MaxQPrefixWidth]uint64
	for limb := 0; limb < rows; limb++ {
		subring := ringQ.SubRings[limb]
		scalarMod := new(big.Int).Mod(scalar, new(big.Int).SetUint64(subring.Modulus)).Uint64()
		scalarMontgomery[limb] = ring.MForm(scalarMod, subring.Modulus, subring.BRedConstant)
	}
	for d := range op0.Value {
		for limb := 0; limb < rows; limb++ {
			subring := ringQ.SubRings[limb]
			subring.MulScalarMontgomery(op0.Value[d].Coeffs[limb], scalarMontgomery[limb], opOut.Value[d].Coeffs[limb])
		}
	}
	*opOut.MetaData = *op0.MetaData
	return nil
}

// Automorphism applies a Fast automorphism to a degree-one ciphertext. It
// operates only on q0 and q1 and does not use evaluation keys or relinearize.
// The input and output may alias. The ciphertexts must have the same level;
// Fast does not perform level reduction as part of this operation.
func (eval *Evaluator) Automorphism(ctIn, ctOut *rlwe.Ciphertext, galEl uint64) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if ctIn == nil || ctOut == nil {
		return errors.New("ctIn and ctOut cannot be nil")
	}
	if ctIn.MetaData == nil || ctOut.MetaData == nil {
		return errors.New("ctIn and ctOut metadata cannot be nil")
	}
	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		return errors.New("Fast automorphism requires degree-one ciphertexts")
	}
	if eval.Parameters.RingType() != ring.Standard {
		return fmt.Errorf("Fast automorphism requires the Standard ring, got %s", eval.Parameters.RingType())
	}
	if ctIn.Level() != ctOut.Level() || ctIn.Level() < 1 {
		return errors.New("Fast automorphism requires equal levels containing q0 and q1")
	}
	if ctIn.N() != eval.Parameters.N() || ctOut.N() != eval.Parameters.N() {
		return errors.New("ciphertext dimensions do not match Fast evaluator parameters")
	}
	if ctIn.IsNTT != ctOut.IsNTT || ctIn.IsMontgomery != ctOut.IsMontgomery {
		return errors.New("Fast automorphism requires equal input/output domains and Montgomery representations")
	}
	for name, ct := range map[string]*rlwe.Ciphertext{"ctIn": ctIn, "ctOut": ctOut} {
		for d := 0; d <= 1; d++ {
			if len(ct.Value[d].Coeffs) < 2 || len(ct.Value[d].Coeffs[0]) != eval.Parameters.N() || len(ct.Value[d].Coeffs[1]) != eval.Parameters.N() {
				return fmt.Errorf("%s has invalid q0/q1 storage", name)
			}
		}
	}

	ringQ := eval.Parameters.RingQ().AtLevel(ctIn.Level())
	rows := maintainedLimbCountForRingAtLevel(ringQ, ctIn.Level())
	if err := validatePrefixRows(ringQ, ctIn.Level(), rows,
		ctIn.Value[0], ctOut.Value[0], ctIn.Value[1], ctOut.Value[1],
	); err != nil {
		return fmt.Errorf("Fast automorphism rows: %w", err)
	}
	if eval.automorphismCore == nil {
		return errors.New("Fast automorphism core is not initialized")
	}
	if err := eval.automorphismCore.ApplyRows(ringQ, galEl, ctIn.IsNTT, rows,
		fastcore.PolynomialPair{Input: ctIn.Value[0], Output: ctOut.Value[0]},
		fastcore.PolynomialPair{Input: ctIn.Value[1], Output: ctOut.Value[1]},
	); err != nil {
		return fmt.Errorf("FastAutomorphism: %w", err)
	}
	*ctOut.MetaData = *ctIn.MetaData
	return nil
}

// AutomorphismQPrefixRows applies a Fast automorphism to exactly rows
// authoritative Q-prefix limbs. It does not use evaluation keys or change
// logical Level/Scale.
func (eval *Evaluator) AutomorphismQPrefixRows(ctIn, ctOut *rlwe.Ciphertext, galEl uint64, rows int) error {
	if eval == nil || ctIn == nil || ctOut == nil || ctIn.MetaData == nil || ctOut.MetaData == nil {
		return errors.New("Fast Q-prefix automorphism evaluator and ciphertext metadata cannot be nil")
	}
	if ctIn.Degree() != 1 || ctOut.Degree() != 1 || ctIn.Level() != ctOut.Level() {
		return errors.New("Fast Q-prefix automorphism requires degree-one ciphertexts at equal levels")
	}
	if eval.Parameters.RingType() != ring.Standard {
		return fmt.Errorf("Fast Q-prefix automorphism requires the Standard ring, got %s", eval.Parameters.RingType())
	}
	if ctIn.N() != eval.Parameters.N() || ctOut.N() != eval.Parameters.N() {
		return errors.New("ciphertext dimensions do not match Fast evaluator parameters")
	}
	if ctIn.IsNTT != ctOut.IsNTT || ctIn.IsMontgomery != ctOut.IsMontgomery {
		return errors.New("Fast Q-prefix automorphism requires matching input/output domains")
	}
	ringQ := eval.Parameters.RingQ().AtLevel(ctIn.Level())
	for name, ct := range map[string]*rlwe.Ciphertext{"ctIn": ctIn, "ctOut": ctOut} {
		for d := 0; d <= 1; d++ {
			if err := validatePrefixRows(ringQ, ctIn.Level(), rows, ct.Value[d]); err != nil {
				return fmt.Errorf("%s component %d: %w", name, d, err)
			}
		}
	}
	if eval.automorphismCore == nil {
		return errors.New("Fast automorphism core is not initialized")
	}
	if err := eval.automorphismCore.ApplyRows(ringQ, galEl, ctIn.IsNTT, rows,
		fastcore.PolynomialPair{Input: ctIn.Value[0], Output: ctOut.Value[0]},
		fastcore.PolynomialPair{Input: ctIn.Value[1], Output: ctOut.Value[1]},
	); err != nil {
		return fmt.Errorf("Fast Q-prefix automorphism: %w", err)
	}
	*ctOut.MetaData = *ctIn.MetaData
	return nil
}

// RotateQPrefixRows applies the slot rotation using exactly rows
// authoritative Q-prefix limbs.
func (eval *Evaluator) RotateQPrefixRows(ctIn, ctOut *rlwe.Ciphertext, k, rows int) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	return eval.AutomorphismQPrefixRows(ctIn, ctOut, eval.Parameters.GaloisElementForRotation(k), rows)
}

// ConjugateQPrefixRows applies complex conjugation using exactly rows
// authoritative Q-prefix limbs.
func (eval *Evaluator) ConjugateQPrefixRows(ctIn, ctOut *rlwe.Ciphertext, rows int) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	return eval.AutomorphismQPrefixRows(ctIn, ctOut, eval.Parameters.GaloisElementForComplexConjugation(), rows)
}

// Rotate applies the Fast automorphism corresponding to a CKKS slot rotation.
func (eval *Evaluator) Rotate(ctIn, ctOut *rlwe.Ciphertext, k int) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	return eval.Automorphism(ctIn, ctOut, eval.Parameters.GaloisElementForRotation(k))
}

// RotateNew returns a newly allocated ciphertext containing the Fast rotation
// of ctIn.
func (eval *Evaluator) RotateNew(ctIn *rlwe.Ciphertext, k int) (*rlwe.Ciphertext, error) {
	if eval == nil || ctIn == nil {
		return nil, errors.New("evaluator and ciphertext cannot be nil")
	}
	ctOut := NewCiphertext(eval.Parameters, 1, ctIn.Level())
	ctOut.IsNTT = ctIn.IsNTT
	ctOut.IsMontgomery = ctIn.IsMontgomery
	return ctOut, eval.Rotate(ctIn, ctOut, k)
}

// Conjugate applies the CKKS complex-conjugation automorphism without an
// evaluation key. Fast supports this only for the Standard ring, through the
// same q0/q1 automorphism path as Rotate.
func (eval *Evaluator) Conjugate(ctIn, ctOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if eval.Parameters.RingType() != ring.Standard {
		return fmt.Errorf("Fast Conjugate requires the Standard ring, got %s", eval.Parameters.RingType())
	}
	return eval.Automorphism(ctIn, ctOut, eval.Parameters.GaloisElementForComplexConjugation())
}

// ConjugateNew returns the Fast conjugate of ctIn.
func (eval *Evaluator) ConjugateNew(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil || ctIn == nil {
		return nil, errors.New("evaluator and ciphertext cannot be nil")
	}
	ctOut := NewCiphertext(eval.Parameters, 1, ctIn.Level())
	ctOut.IsNTT = ctIn.IsNTT
	ctOut.IsMontgomery = ctIn.IsMontgomery
	return ctOut, eval.Conjugate(ctIn, ctOut)
}
