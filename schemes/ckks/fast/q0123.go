package fast

import (
	"math/bits"
)

// mul192By64 multiplies a fixed-width value by a uint64 and reports whether
// the mathematical product exceeds the 192-bit representation.
func mul192By64(value uint192, scalar uint64) (product uint192, overflow bool) {
	hi0, lo0 := bits.Mul64(value.lo, scalar)
	hi1, lo1 := bits.Mul64(value.mid, scalar)
	hi2, lo2 := bits.Mul64(value.hi, scalar)

	mid, carry := bits.Add64(hi0, lo1, 0)
	hi, carryOut := bits.Add64(hi1, lo2, carry)
	return uint192{lo: lo0, mid: mid, hi: hi}, hi2 != 0 || carryOut != 0
}

func add192(left, right uint192) (sum uint192, overflow bool) {
	lo, carry := bits.Add64(left.lo, right.lo, 0)
	mid, carry := bits.Add64(left.mid, right.mid, carry)
	hi, carry := bits.Add64(left.hi, right.hi, carry)
	return uint192{lo: lo, mid: mid, hi: hi}, carry != 0
}

func half192(value uint192) uint192 {
	return uint192{
		lo:  (value.lo >> 1) | (value.mid << 63),
		mid: (value.mid >> 1) | (value.hi << 63),
		hi:  value.hi >> 1,
	}
}

// q0123Modulus computes the exact product q0*q1*q2*q3 and floor(product/2).
func q0123Modulus(q0, q1, q2, q3 uint64) (modulus, half uint192, q012 uint192, overflow bool) {
	_, _, q01Lo, q01Hi := q012Modulus(q0, q1, q2)
	q012 = uint192{lo: q01Lo, mid: q01Hi}
	q012, overflow = mul192By64(q012, q2)
	if overflow {
		return uint192{}, uint192{}, uint192{}, true
	}
	modulus, overflow = mul192By64(q012, q3)
	if overflow {
		return uint192{}, uint192{}, q012, true
	}
	return modulus, half192(modulus), q012, false
}

// crtQ0123 reconstructs the canonical residue in [0,q0*q1*q2*q3) using
// Garner reconstruction. The caller validates the moduli, inverse, and
// supported 192-bit product range before entering a coefficient loop.
func crtQ0123(r0, r1, r2, r3, q0, q1, q2, q3, q0InvQ1, q01InvQ2, q012InvQ3 uint64) uint192 {
	x012 := crtQ012(r0, r1, r2, q0, q1, q2, q0InvQ1, q01InvQ2)
	x012ModQ3 := mod192By64(x012, q3)
	var delta uint64
	if r3 >= x012ModQ3 {
		delta = r3 - x012ModQ3
	} else {
		delta = q3 - (x012ModQ3 - r3)
	}
	prodHi, prodLo := bits.Mul64(delta, q012InvQ3)
	_, t := bits.Div64(prodHi, prodLo, q3)

	_, _, q01Lo, q01Hi := q012Modulus(q0, q1, q2)
	q012, _ := mul192By64(uint192{lo: q01Lo, mid: q01Hi}, q2)
	term, _ := mul192By64(q012, t)
	x, _ := add192(x012, term)
	return x
}
