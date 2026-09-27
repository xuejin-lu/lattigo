package fast

import (
	"errors"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// MulThenAddOneBitScalarGuard performs one explicitly requested guarded scalar
// accumulation. The accumulator is promoted by one binary bit, the ordinary
// MulThenAdd implementation is used unchanged, and the result is contracted
// with centered rounded division by two. This is deliberately not part of the
// generic MulThenAdd policy.
func (eval *Evaluator) MulThenAddOneBitScalarGuard(op0 *rlwe.Ciphertext, op1 rlwe.Operand, opOut *rlwe.Ciphertext) error {
	if eval == nil || op0 == nil || opOut == nil {
		return errors.New("Fast one-bit scalar guard evaluator and operands cannot be nil")
	}
	level := opOut.Level()
	if op0.Level() < level {
		level = op0.Level()
	}
	return eval.MulThenAddOneBitScalarGuardQPrefixRows(op0, op1, maintainedLimbCount(&eval.Parameters, level), opOut)
}

// MulThenAddOneBitScalarGuardQPrefixRows performs the explicit guarded scalar
// accumulation and centered contraction over exactly rows authoritative Q
// prefix residues.
func (eval *Evaluator) MulThenAddOneBitScalarGuardQPrefixRows(op0 *rlwe.Ciphertext, op1 rlwe.Operand, rows int, opOut *rlwe.Ciphertext) error {
	if eval == nil || op0 == nil || opOut == nil {
		return errors.New("Fast one-bit scalar guard evaluator and operands cannot be nil")
	}
	if opOut.MetaData == nil || op0.MetaData == nil {
		return errors.New("Fast one-bit scalar guard metadata cannot be nil")
	}
	if opOut.Level() < 1 || op0.Level() < 1 {
		return errors.New("Fast one-bit scalar guard requires q0 and q1")
	}
	if !opOut.IsNTT || !op0.IsNTT || opOut.IsMontgomery != op0.IsMontgomery {
		return errors.New("Fast one-bit scalar guard requires matching NTT/Montgomery operands")
	}
	if opOut.Degree() != op0.Degree() {
		return errors.New("Fast one-bit scalar guard requires matching degrees")
	}
	level := opOut.Level()
	if op0.Level() < level {
		level = op0.Level()
	}
	if err := eval.validateExplicitRows(level, rows); err != nil {
		return err
	}
	for d := range opOut.Value {
		if err := validatePrefixRows(eval.Parameters.RingQ(), opOut.Level(), rows, opOut.Value[d]); err != nil {
			return err
		}
	}
	for d := range op0.Value {
		if err := validatePrefixRows(eval.Parameters.RingQ(), op0.Level(), rows, op0.Value[d]); err != nil {
			return err
		}
	}

	nativeLevel := opOut.Level()
	nativeDegree := opOut.Degree()
	nativeNTT := opOut.IsNTT
	nativeMontgomery := opOut.IsMontgomery
	nativeMeta := *opOut.MetaData

	if err := eval.MulIntegerQPrefixRows(opOut, big.NewInt(2), rows, opOut); err != nil {
		return fmtGuardError("promote accumulator", err)
	}
	nativeScale := opOut.Scale
	opOut.Scale = nativeScale.Mul(rlwe.NewScale(2))
	guardedScale := opOut.Scale

	if err := eval.MulThenAddQPrefixRows(op0, op1, rows, opOut); err != nil {
		return fmtGuardError("scalar MulThenAdd", err)
	}
	if opOut.Level() != nativeLevel || opOut.Degree() != nativeDegree || opOut.IsNTT != nativeNTT || opOut.IsMontgomery != nativeMontgomery {
		return errors.New("Fast one-bit scalar guard changed ciphertext structure before contraction")
	}

	if err := eval.contractCenteredRoundedDivideByTwoRows(opOut, rows); err != nil {
		return fmtGuardError("centered rounded divide by two", err)
	}
	if opOut.Level() != nativeLevel || opOut.Degree() != nativeDegree || opOut.IsNTT != nativeNTT || opOut.IsMontgomery != nativeMontgomery {
		return errors.New("Fast one-bit scalar guard changed ciphertext structure during contraction")
	}
	*opOut.MetaData = nativeMeta
	opOut.Scale = guardedScale.Div(rlwe.NewScale(2))
	if opOut.Level() != nativeLevel || opOut.Degree() != nativeDegree || !opOut.Scale.Equal(guardedScale.Div(rlwe.NewScale(2))) {
		return errors.New("Fast one-bit scalar guard failed to restore native metadata")
	}
	return nil
}

// fmtGuardError keeps the production errors contextual without importing the
// formatting package into the fixed-width arithmetic implementation.
func fmtGuardError(operation string, err error) error {
	return errors.New("Fast one-bit scalar guard " + operation + ": " + err.Error())
}

func (eval *Evaluator) contractCenteredRoundedDivideByTwo(ct *rlwe.Ciphertext) error {
	if eval == nil || ct == nil {
		return errors.New("Fast centered divide-by-two operands cannot be nil")
	}
	return eval.contractCenteredRoundedDivideByTwoRows(ct, maintainedLimbCount(&eval.Parameters, ct.Level()))
}

func (eval *Evaluator) contractCenteredRoundedDivideByTwoRows(ct *rlwe.Ciphertext, rows int) error {
	if eval == nil || ct == nil || ct.MetaData == nil {
		return errors.New("Fast centered divide-by-two operands cannot be nil")
	}
	if !ct.IsNTT || len(ct.Value) == 0 || ct.Level() < 1 {
		return errors.New("Fast centered divide-by-two requires NTT q0/q1 ciphertexts")
	}
	if err := eval.validateExplicitRows(ct.Level(), rows); err != nil {
		return err
	}

	ringQ := eval.Parameters.RingQ()
	scratch := &eval.rescaleScratch
	for component := range ct.Value {
		if err := prefixToCoefficientRows(ringQ, ct.Value[component], rows, true, ct.IsMontgomery, scratch.coeff); err != nil {
			return err
		}
		for k := 0; k < ringQ.N(); k++ {
			var residues [MaxQPrefixWidth]uint64
			for row := 0; row < rows; row++ {
				residues[row] = scratch.coeff.Coeffs[row][k]
			}
			value := reconstructQPrefix(rows, residues, scratch)
			magnitude, negative := centeredQPrefix(value, scratch.modulus[rows-1], scratch.half[rows-1])
			magnitude = roundedMagnitude192ByTwo(magnitude)
			for row := 0; row < rows; row++ {
				scratch.result.Coeffs[row][k] = signedResidue192(magnitude, negative, scratch.q[row])
			}
		}
		for row := 0; row < rows; row++ {
			ringQ.SubRings[row].NTT(scratch.result.Coeffs[row], ct.Value[component].Coeffs[row])
			if ct.IsMontgomery {
				ringQ.SubRings[row].MForm(ct.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
			}
		}
	}
	return nil
}

func roundedMagnitude192ByTwo(value uint192) uint192 {
	quotient := uint192{
		lo:  (value.lo >> 1) | (value.mid << 63),
		mid: (value.mid >> 1) | (value.hi << 63),
		hi:  value.hi >> 1,
	}
	if value.lo&1 != 0 {
		var carry uint64
		quotient.lo, carry = bits.Add64(quotient.lo, 1, 0)
		quotient.mid, carry = bits.Add64(quotient.mid, 0, carry)
		quotient.hi, _ = bits.Add64(quotient.hi, 0, carry)
	}
	return quotient
}

func roundedMagnitude128ByTwo(xLo, xHi, qLo, qHi, halfLo, halfHi uint64) (uint64, uint64, bool) {
	negative := xHi > halfHi || (xHi == halfHi && xLo > halfLo)
	if negative {
		borrow := uint64(0)
		qLo, borrow = bits.Sub64(qLo, xLo, 0)
		qHi, _ = bits.Sub64(qHi, xHi, borrow)
	} else {
		qLo, qHi = xLo, xHi
	}

	// qLo/qHi is the non-negative centered magnitude. Symmetric nearest
	// integer division rounds every odd magnitude away from zero.
	resultLo := qLo >> 1
	resultHi := qHi >> 1
	resultLo |= qHi << 63
	if qLo&1 != 0 {
		var carry uint64
		resultLo, carry = bits.Add64(resultLo, 1, 0)
		resultHi += carry
	}
	return resultLo, resultHi, negative
}

func signedResidue128(magnitudeLo, magnitudeHi uint64, negative bool, modulus uint64) uint64 {
	_, remainder := bits.Div64(0, magnitudeHi, modulus)
	_, remainder = bits.Div64(remainder, magnitudeLo, modulus)
	if negative && remainder != 0 {
		return modulus - remainder
	}
	return remainder
}
