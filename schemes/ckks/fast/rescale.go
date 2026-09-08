package fast

import (
	"errors"
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const (
	fastRescaleMaxQ0Bits      = 55
	fastRescaleMaxQ1Bits      = 39
	fastRescaleMaxQ01Bits     = 94
	fastRescaleMinDivisorBits = 32
)

// Rescale applies the Standard CKKS rescale semantics using only the
// actively-maintained q0/q1 residues. The input must be NTT-domain and
// Montgomery-form, as produced by the Standard CKKS evaluator.
func (eval *Evaluator) Rescale(op0, opOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if op0 == nil || opOut == nil {
		return errors.New("op0 and opOut cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return errors.New("op0 and opOut metadata cannot be nil")
	}
	nbRescales := eval.Parameters.LevelsConsumedPerRescaling()
	if op0.Level() < nbRescales {
		return errors.New("cannot Rescale: input Ciphertext level is too low")
	}
	return eval.rescaleN(op0, nbRescales, opOut)
}

// RescaleTo repeatedly rescales until the scale reaches minScale, or another
// rescale would take it below minScale/2. It preserves the Standard stopping
// rule while allowing the output to reach Level 0.
func (eval *Evaluator) RescaleTo(op0 *rlwe.Ciphertext, minScale rlwe.Scale, opOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if op0 == nil || opOut == nil {
		return errors.New("op0 and opOut cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return errors.New("op0 and opOut metadata cannot be nil")
	}
	if minScale.Cmp(rlwe.NewScale(0)) != 1 {
		return errors.New("cannot RescaleTo: minScale is not positive")
	}
	if op0.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return errors.New("cannot RescaleTo: ciphertext scale is not positive")
	}
	if op0.Level() == 0 {
		return errors.New("cannot RescaleTo: input Ciphertext already at level 0")
	}

	threshold := minScale.Div(rlwe.NewScale(2))
	scale := op0.Scale
	newLevel := op0.Level()
	nbRescales := 0
	for newLevel > 0 {
		candidate := scale.Div(rlwe.NewScale(eval.Parameters.Q()[newLevel]))
		if candidate.Cmp(threshold) == -1 {
			break
		}
		scale = candidate
		newLevel--
		nbRescales++
	}

	if nbRescales == 0 {
		if op0 != opOut {
			opOut.Copy(op0)
		}
		return nil
	}
	return eval.rescaleN(op0, nbRescales, opOut)
}

func (eval *Evaluator) rescaleN(op0 *rlwe.Ciphertext, nbRescales int, opOut *rlwe.Ciphertext) error {
	if nbRescales <= 0 || op0.Level() < nbRescales {
		return errors.New("invalid number of Fast rescale levels")
	}
	if op0.N() != eval.Parameters.N() || opOut.N() != eval.Parameters.N() {
		return errors.New("ciphertext dimensions do not match Fast evaluator parameters")
	}
	validate := validateFastRescaleDomain(op0)
	if validate != nil {
		return validate
	}

	current := op0.CopyNew()
	for i := 0; i < nbRescales; i++ {
		level := current.Level()
		var next *rlwe.Ciphertext
		if i == nbRescales-1 {
			next = opOut
			next.Resize(current.Degree(), level-1)
		} else {
			next = ckks.NewCiphertext(eval.Parameters, current.Degree(), level-1)
		}
		if err := fastRescaleOnce(eval.Parameters.RingQ().AtLevel(level), current, next); err != nil {
			return fmt.Errorf("Fast rescale at level %d: %w", level, err)
		}
		current = next
	}
	return nil
}

func fastRescaleOnce(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext) error {
	if op0.Level() < 1 {
		return errors.New("Fast Rescale requires input Level >= 1")
	}
	if err := validateFastRescaleDomain(op0); err != nil {
		return err
	}
	outMontgomery := op0.IsMontgomery
	if opOut.Degree() != op0.Degree() || opOut.Level() != op0.Level()-1 {
		return errors.New("Fast Rescale output must have the same degree and one lower level")
	}

	q0 := ringQ.SubRings[0].Modulus
	q1 := ringQ.SubRings[1].Modulus
	d := ringQ.SubRings[ringQ.Level()].Modulus
	inv, ok := inverseMod(q0%q1, q1)
	if !ok {
		return errors.New("Fast Rescale requires coprime q0 and q1")
	}
	if err := validateFastRescaleRange(q0, q1, d); err != nil {
		return err
	}

	q01Hi, q01Lo := bits.Mul64(q0, q1)
	halfLo := (q01Lo >> 1) | (q01Hi << 63)
	halfHi := q01Hi >> 1
	results := make([]ring.Poly, op0.Degree()+1)
	for component := range op0.Value {
		coeff := ring.NewPoly(ringQ.N(), 1)
		if err := FastPartialINTT(ringQ, op0.Value[component], coeff); err != nil {
			return err
		}
		if op0.IsMontgomery {
			for limb := 0; limb < 2; limb++ {
				ringQ.SubRings[limb].IMForm(coeff.Coeffs[limb], coeff.Coeffs[limb])
			}
		}
		result := ring.NewPoly(ringQ.N(), opOut.Level())
		for k := 0; k < ringQ.N(); k++ {
			xLo, xHi := crtQ01(coeff.Coeffs[0][k], coeff.Coeffs[1][k], q0, q1, inv)
			y, negative := roundedMagnitude128(xLo, xHi, q01Lo, q01Hi, halfLo, halfHi, d)
			result.Coeffs[0][k] = signedResidue(y, negative, q0)
			if opOut.Level() >= 1 {
				result.Coeffs[1][k] = signedResidue(y, negative, q1)
			}
		}

		if opOut.Level() >= 1 {
			nttResult := ring.NewPoly(ringQ.N(), opOut.Level())
			if err := FastPartialNTT(ringQ, result, nttResult); err != nil {
				return err
			}
			result = nttResult
		} else {
			nttResult := ring.NewPoly(ringQ.N(), 0)
			ringQ.SubRings[0].NTT(result.Coeffs[0], nttResult.Coeffs[0])
			result = nttResult
		}
		if outMontgomery {
			for limb := 0; limb <= opOut.Level() && limb < 2; limb++ {
				ringQ.SubRings[limb].MForm(result.Coeffs[limb], result.Coeffs[limb])
			}
		}
		results[component] = result
	}

	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale.Div(rlwe.NewScale(d))
	opOut.Resize(op0.Degree(), op0.Level()-1)
	for component := range results {
		opOut.Value[component] = results[component]
	}
	return nil
}

func validateFastRescaleDomain(ct *rlwe.Ciphertext) error {
	if !ct.IsNTT {
		return errors.New("Fast Rescale requires NTT-domain ciphertexts")
	}
	if len(ct.Value) == 0 || ct.Value[0].Level() < 1 {
		return errors.New("Fast Rescale requires maintained q0 and q1 storage")
	}
	return nil
}

func validateFastRescaleRange(q0, q1, divisor uint64) error {
	if bits.Len64(q0) > fastRescaleMaxQ0Bits || bits.Len64(q1) > fastRescaleMaxQ1Bits {
		return fmt.Errorf("unsupported Fast Rescale moduli: q0 must be <= %d bits and q1 <= %d bits", fastRescaleMaxQ0Bits, fastRescaleMaxQ1Bits)
	}
	q01Hi, q01Lo := bits.Mul64(q0, q1)
	if q01Hi != 0 && bits.Len64(q01Hi)+64 > fastRescaleMaxQ01Bits || q01Hi == 0 && bits.Len64(q01Lo) > fastRescaleMaxQ01Bits {
		return fmt.Errorf("unsupported Fast Rescale q0*q1 range: requires product < 2^%d", fastRescaleMaxQ01Bits)
	}
	if bits.Len64(divisor) < fastRescaleMinDivisorBits {
		return fmt.Errorf("unsupported Fast Rescale divisor: requires at least %d bits", fastRescaleMinDivisorBits)
	}
	return nil
}

func inverseMod(a, modulus uint64) (uint64, bool) {
	if modulus == 0 || a == 0 {
		return 0, false
	}
	var oldR, r = int64(a), int64(modulus)
	var oldT, t int64 = 1, 0
	for r != 0 {
		q := oldR / r
		oldR, r = r, oldR-q*r
		oldT, t = t, oldT-q*t
	}
	if oldR != 1 {
		return 0, false
	}
	if oldT < 0 {
		oldT += int64(modulus)
	}
	return uint64(oldT), true
}

func crtQ01(r0, r1, q0, q1, inverse uint64) (lo, hi uint64) {
	r0mod := r0 % q1
	var delta uint64
	if r1 >= r0mod {
		delta = r1 - r0mod
	} else {
		delta = q1 - (r0mod - r1)
	}
	prodHi, prodLo := bits.Mul64(delta, inverse)
	_, t := bits.Div64(prodHi, prodLo, q1)
	prodHi, prodLo = bits.Mul64(t, q0)
	lo, carry := bits.Add64(prodLo, r0, 0)
	hi, _ = bits.Add64(prodHi, 0, carry)
	return lo, hi
}

func roundedMagnitude128(xLo, xHi, qLo, qHi, halfLo, halfHi, divisor uint64) (uint64, bool) {
	negative := xHi > halfHi || (xHi == halfHi && xLo >= halfLo)
	if negative {
		borrow := uint64(0)
		qLo, borrow = bits.Sub64(qLo, xLo, 0)
		qHi, _ = bits.Sub64(qHi, xHi, borrow)
	} else {
		qLo, qHi = xLo, xHi
	}
	_, remainder := bits.Div64(0, qHi, divisor)
	quotient, remainder := bits.Div64(remainder, qLo, divisor)
	if remainder > divisor/2 {
		quotient++
	}
	return quotient, negative
}

func signedResidue(magnitude uint64, negative bool, modulus uint64) uint64 {
	r := magnitude % modulus
	if negative && r != 0 {
		return modulus - r
	}
	return r
}
