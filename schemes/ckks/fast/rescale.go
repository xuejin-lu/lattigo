package fast

import (
	"errors"
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

const (
	// The LogN=16 public parameter builder can generate a 56-bit q0 and a
	// 39-bit q1. For this supported domain q0 < 2^56 and q1 < 2^39, hence
	// q0*q1 < 2^95 and the two-word CRT product remains below 128 bits.
	// The centered magnitude is < 2^94. Since the divisor validator requires
	// d >= 2^31, roundedMagnitude128 always satisfies hi < d for its first
	// bits.Div64 call; the second call receives a remainder < d by contract.
	// In crtQ01, delta*inverse < q1^2 < 2^78, so its high word is below the
	// actual 39-bit q1 used by the supported profile.
	fastRescaleMaxQ0Bits      = 56
	fastRescaleMaxQ1Bits      = 39
	fastRescaleMaxQ01Bits     = 95
	fastRescaleMinDivisorBits = 32
)

type fastRescaleScratch struct {
	coeff, result ring.Poly
	staged        []ring.Poly
	diagParent    uint64
	q             [MaxQPrefixWidth]uint64
	inverse       [MaxQPrefixWidth]uint64
	modulus       [MaxQPrefixWidth]uint192
	half          [MaxQPrefixWidth]uint192
	initErr       error
	widthErr      [MaxQPrefixWidth]error
	// Scalar guards still use the legacy Q01 CRT scratch; Rescale itself uses
	// the explicit prefix arrays above.
	q0, q1, q0InverseModQ1 uint64
	q01Lo, q01Hi           uint64
	halfLo, halfHi         uint64
}

func newFastRescaleScratch(ringQ *ring.Ring) fastRescaleScratch {
	scratch := fastRescaleScratch{}
	if ringQ == nil || ringQ.Level() < 1 {
		return scratch
	}
	prefixWidth := qPrefixWidthOrPanic(ringQ.Level())
	scratch.coeff = ring.NewPoly(ringQ.N(), prefixWidth-1)
	scratch.result = ring.NewPoly(ringQ.N(), prefixWidth-1)
	// A degree-one ciphertext is the common Rescale case. Allocate its staged
	// component buffers with the evaluator so steady-state Rescale does not
	// allocate coefficient backing on every call. Higher degrees grow lazily.
	scratch.staged = make([]ring.Poly, 2)
	for i := range scratch.staged {
		scratch.staged[i] = ring.NewPoly(ringQ.N(), prefixWidth-1)
	}
	if len(ringQ.SubRings) < prefixWidth {
		scratch.initErr = fmt.Errorf("Fast Rescale requires %d configured q subrings", prefixWidth)
		return scratch
	}
	product := uint192{lo: 1}
	productValid := true
	for row := 0; row < prefixWidth; row++ {
		q := ringQ.SubRings[row].Modulus
		if q < 3 || q&1 == 0 || q > uint64(^uint64(0)>>1) {
			scratch.widthErr[row] = fmt.Errorf("unsupported Fast Rescale modulus q%d=%d", row, q)
			productValid = false
			continue
		}
		scratch.q[row] = q
		if !productValid {
			scratch.widthErr[row] = fmt.Errorf("Fast Rescale Q-prefix through q%d is unavailable", row)
			continue
		}
		var overflow bool
		product, overflow = mul192By64(product, q)
		if overflow {
			scratch.widthErr[row] = fmt.Errorf("Fast Rescale Q-prefix product through q%d exceeds 192 bits", row)
			productValid = false
			continue
		}
		scratch.modulus[row] = product
		scratch.half[row] = half192(product)
		if row > 0 {
			inverse, ok := inverseMod(mod192By64(scratch.modulus[row-1], q), q)
			if !ok {
				scratch.widthErr[row] = fmt.Errorf("Fast Rescale Q-prefix through q%d is not pairwise coprime", row)
				productValid = false
				continue
			}
			scratch.inverse[row] = inverse
		}
	}
	if prefixWidth >= 2 {
		scratch.q0, scratch.q1 = scratch.q[0], scratch.q[1]
		scratch.q0InverseModQ1 = scratch.inverse[1]
		scratch.q01Hi, scratch.q01Lo = bits.Mul64(scratch.q0, scratch.q1)
		scratch.halfLo = (scratch.q01Lo >> 1) | (scratch.q01Hi << 63)
		scratch.halfHi = scratch.q01Hi >> 1
	}
	return scratch
}

func (scratch *fastRescaleScratch) validateWidth(rows int) error {
	if scratch == nil {
		return errors.New("Fast Rescale scratch cannot be nil")
	}
	if scratch.initErr != nil {
		return scratch.initErr
	}
	if rows < 1 || rows > MaxQPrefixWidth {
		return fmt.Errorf("invalid Fast Rescale Q-prefix width %d", rows)
	}
	return scratch.widthErr[rows-1]
}

func (scratch *fastRescaleScratch) ensureStagedComponents(ringQ *ring.Ring, components int) error {
	if scratch == nil || ringQ == nil {
		return errors.New("Fast Rescale staging requires non-nil scratch and ring")
	}
	if components < 1 {
		return errors.New("Fast Rescale staging requires at least one component")
	}
	width := len(scratch.coeff.Coeffs)
	for len(scratch.staged) < components {
		scratch.staged = append(scratch.staged, ring.NewPoly(ringQ.N(), width-1))
	}
	for component := range scratch.staged {
		if len(scratch.staged[component].Coeffs) < width || scratch.staged[component].N() != ringQ.N() {
			scratch.staged[component] = ring.NewPoly(ringQ.N(), width-1)
		}
	}
	return nil
}

// Rescale applies Standard CKKS rescale semantics using the fixed Q-prefix
// policy w_Q(Level)=min(Level+1, 4). The logical top modulus is the divisor
// even when it is outside that prefix. The input must be NTT-domain.
// Both ordinary and Montgomery-form NTT representations are supported and the
// output preserves the input representation.
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
	sourceRows, err := QPrefixWidth(op0.Level())
	if err != nil {
		return err
	}
	return eval.rescaleNQPrefix(op0, nbRescales, sourceRows, opOut)
}

// RescaleQPrefixRows applies Standard CKKS Rescale semantics using exactly
// rows authoritative source limbs. Rows outside the resulting source-derived
// target prefix remain unavailable in the output backing.
func (eval *Evaluator) RescaleQPrefixRows(op0 *rlwe.Ciphertext, rows int, opOut *rlwe.Ciphertext) error {
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
	width, err := QPrefixWidth(op0.Level())
	if err != nil {
		return err
	}
	if rows < 1 || rows > width {
		return fmt.Errorf("Fast Rescale row count %d must be in [1,%d] at level %d", rows, width, op0.Level())
	}
	return eval.rescaleNQPrefix(op0, nbRescales, rows, opOut)
}

// RescaleTo repeatedly rescales until the scale reaches minScale, or another
// rescale would take it below minScale/2. It preserves the Standard stopping
// rule while allowing the output to reach Level 0.
func (eval *Evaluator) RescaleTo(op0 *rlwe.Ciphertext, minScale rlwe.Scale, opOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if op0 == nil {
		return errors.New("op0 cannot be nil")
	}
	sourceRows, err := QPrefixWidth(op0.Level())
	if err != nil {
		return err
	}
	return eval.rescaleToNQPrefix(op0, minScale, sourceRows, opOut)
}

// RescaleToQPrefixRows applies the Standard RescaleTo stopping rule using an
// explicit source-prefix authority. A no-op preserves Level and Scale while
// retaining only the explicitly authorized rows.
func (eval *Evaluator) RescaleToQPrefixRows(op0 *rlwe.Ciphertext, minScale rlwe.Scale, rows int, opOut *rlwe.Ciphertext) error {
	if eval == nil || op0 == nil || opOut == nil {
		return errors.New("Fast RescaleTo evaluator and ciphertexts cannot be nil")
	}
	if op0.Level() < 1 {
		return errors.New("cannot RescaleTo: input Ciphertext already at level 0")
	}
	if err := eval.validateExplicitRows(op0.Level(), rows); err != nil {
		return err
	}
	return eval.rescaleToNQPrefix(op0, minScale, rows, opOut)
}

// rescaleToNQPrefix applies RescaleTo's Standard stopping rule while keeping
// source authority explicit. Production wrappers pass the fixed Q-prefix
// width; explicit row-authority callers may pass a narrower supported prefix.
func (eval *Evaluator) rescaleToNQPrefix(op0 *rlwe.Ciphertext, minScale rlwe.Scale, sourceRows int, opOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if op0 == nil || opOut == nil {
		return errors.New("op0 and opOut cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return errors.New("op0 and opOut metadata cannot be nil")
	}
	if eval.rescaleCore == nil {
		return errors.New("Fast Rescale core is not initialized")
	}
	return eval.rescaleCore.ApplyToRows(eval.Parameters.RingQ(), op0, minScale, sourceRows, opOut)
}

func validateFastRescaleDomain(ct *rlwe.Ciphertext) error {
	if !ct.IsNTT {
		return errors.New("Fast Rescale requires NTT-domain ciphertexts")
	}
	if len(ct.Value) == 0 || ct.Value[0].Level() < 1 {
		return errors.New("Fast Rescale requires a non-zero logical Level")
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
	// For odd Q01, floor(Q01/2) is still the positive centered endpoint.
	negative := xHi > halfHi || (xHi == halfHi && xLo > halfLo)
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

func roundedMagnitude64(magnitude, divisor uint64) uint64 {
	quotient, remainder := magnitude/divisor, magnitude%divisor
	if remainder > divisor/2 {
		quotient++
	}
	return quotient
}

func signedResidue(magnitude uint64, negative bool, modulus uint64) uint64 {
	r := magnitude % modulus
	if negative && r != 0 {
		return modulus - r
	}
	return r
}
