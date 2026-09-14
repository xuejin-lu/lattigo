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

	nativeLevel := opOut.Level()
	nativeDegree := opOut.Degree()
	nativeNTT := opOut.IsNTT
	nativeMontgomery := opOut.IsMontgomery
	nativeMeta := *opOut.MetaData

	if err := eval.MulIntegerMaintained(opOut, big.NewInt(2), opOut); err != nil {
		return fmtGuardError("promote accumulator", err)
	}
	nativeScale := opOut.Scale
	opOut.Scale = nativeScale.Mul(rlwe.NewScale(2))
	guardedScale := opOut.Scale

	if err := eval.MulThenAdd(op0, op1, opOut); err != nil {
		return fmtGuardError("scalar MulThenAdd", err)
	}
	if opOut.Level() != nativeLevel || opOut.Degree() != nativeDegree || opOut.IsNTT != nativeNTT || opOut.IsMontgomery != nativeMontgomery {
		return errors.New("Fast one-bit scalar guard changed ciphertext structure before contraction")
	}

	if err := eval.contractCenteredRoundedDivideByTwo(opOut); err != nil {
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
	if eval == nil || ct == nil || ct.MetaData == nil {
		return errors.New("Fast centered divide-by-two operands cannot be nil")
	}
	if !ct.IsNTT || len(ct.Value) == 0 || ct.Level() < 1 {
		return errors.New("Fast centered divide-by-two requires NTT q0/q1 ciphertexts")
	}
	if err := validateFastRescaleRange(eval.Parameters.RingQ().SubRings[0].Modulus, eval.Parameters.RingQ().SubRings[1].Modulus, 1<<32); err != nil {
		return err
	}

	ringQ := eval.Parameters.RingQ()
	scratch := &eval.rescaleScratch
	for component := range ct.Value {
		if err := FastPartialINTT(ringQ, ct.Value[component], scratch.coeff); err != nil {
			return err
		}
		if ct.IsMontgomery {
			for limb := 0; limb < 2; limb++ {
				ringQ.SubRings[limb].IMForm(scratch.coeff.Coeffs[limb], scratch.coeff.Coeffs[limb])
			}
		}
		for k := 0; k < ringQ.N(); k++ {
			xLo, xHi := crtQ01(scratch.coeff.Coeffs[0][k], scratch.coeff.Coeffs[1][k], scratch.q0, scratch.q1, scratch.q0InverseModQ1)
			magnitudeLo, magnitudeHi, negative := roundedMagnitude128ByTwo(xLo, xHi, scratch.q01Lo, scratch.q01Hi, scratch.halfLo, scratch.halfHi)
			scratch.result.Coeffs[0][k] = signedResidue128(magnitudeLo, magnitudeHi, negative, scratch.q0)
			scratch.result.Coeffs[1][k] = signedResidue128(magnitudeLo, magnitudeHi, negative, scratch.q1)
		}
		ringQ.SubRings[0].NTT(scratch.result.Coeffs[0], ct.Value[component].Coeffs[0])
		ringQ.SubRings[1].NTT(scratch.result.Coeffs[1], ct.Value[component].Coeffs[1])
		if ct.IsMontgomery {
			for limb := 0; limb < 2; limb++ {
				ringQ.SubRings[limb].MForm(ct.Value[component].Coeffs[limb], ct.Value[component].Coeffs[limb])
			}
		}
	}
	return nil
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
