package fast

import (
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// q012Enabled selects the bounded LogN13 production profile. The profile is
// identified by its generated modulus shape rather than by an application
// flag, so Standard and Fast continue to consume the same frontend config.
func q012Enabled(params rlwe.ParameterProvider) bool {
	if params == nil || params.GetRLWEParameters() == nil {
		return false
	}
	p := params.GetRLWEParameters()
	q := p.Q()
	return len(q) >= 3 && bits.Len64(q[0]) == 56 && bits.Len64(q[1]) <= 39 && bits.Len64(q[2]) <= 40
}

// maintainedLimbCount returns the actively maintained residue count for the
// bounded Fast mode at the requested logical level.
func maintainedLimbCount(params rlwe.ParameterProvider, level int) int {
	count := 2
	if q012Enabled(params) {
		count = 3
	}
	if level+1 < count {
		return level + 1
	}
	return count
}

// MaintainedLimbCount exposes the active residue count to the bounded CKKS
// circuit adapters without exposing any backend-specific mode selector.
func MaintainedLimbCount(params rlwe.ParameterProvider, level int) int {
	return maintainedLimbCount(params, level)
}

func maintainedLimbCountForRing(ringQ *ring.Ring) int {
	if ringQ == nil {
		return 0
	}
	return maintainedLimbCountForRingAtLevel(ringQ, ringQ.Level())
}

func maintainedLimbCountForRingAtLevel(ringQ *ring.Ring, level int) int {
	if ringQ == nil {
		return 0
	}
	if level > ringQ.Level() {
		level = ringQ.Level()
	}
	count := 2
	if level >= 2 && len(ringQ.SubRings) >= 3 && bits.Len64(ringQ.SubRings[0].Modulus) == 56 && bits.Len64(ringQ.SubRings[1].Modulus) <= 39 && bits.Len64(ringQ.SubRings[2].Modulus) <= 40 {
		count = 3
	}
	if ringQ.Level()+1 < count {
		return ringQ.Level() + 1
	}
	return count
}

// uint192 is the fixed-width production representation used by the bounded
// q0/q1/q2 reconstruction path. The supported LogN13 profile is below 135
// bits, so no arbitrary-precision allocation is needed per coefficient.
type uint192 struct {
	lo, mid, hi uint64
}

func add128To192(lo, hi uint64, value uint192) uint192 {
	lo, carry := bits.Add64(value.lo, lo, 0)
	mid, carry := bits.Add64(value.mid, hi, carry)
	hi, _ = bits.Add64(value.hi, 0, carry)
	return uint192{lo: lo, mid: mid, hi: hi}
}

func mul128By64(lo, hi, scalar uint64) uint192 {
	hi0, lo0 := bits.Mul64(lo, scalar)
	hi1, lo1 := bits.Mul64(hi, scalar)
	mid, carry := bits.Add64(hi0, lo1, 0)
	hi2, _ := bits.Add64(hi1, 0, carry)
	return uint192{lo: lo0, mid: mid, hi: hi2}
}

func sub192(left, right uint192) uint192 {
	lo, borrow := bits.Sub64(left.lo, right.lo, 0)
	mid, borrow := bits.Sub64(left.mid, right.mid, borrow)
	hi, _ := bits.Sub64(left.hi, right.hi, borrow)
	return uint192{lo: lo, mid: mid, hi: hi}
}

func cmp192(left, right uint192) int {
	if left.hi != right.hi {
		if left.hi < right.hi {
			return -1
		}
		return 1
	}
	if left.mid != right.mid {
		if left.mid < right.mid {
			return -1
		}
		return 1
	}
	if left.lo < right.lo {
		return -1
	}
	if left.lo > right.lo {
		return 1
	}
	return 0
}

func q012Modulus(q0, q1, q2 uint64) (modulus, half uint192, q01Lo, q01Hi uint64) {
	q01Hi, q01Lo = bits.Mul64(q0, q1)
	modulus = mul128By64(q01Lo, q01Hi, q2)
	half = uint192{lo: (modulus.lo >> 1) | (modulus.mid << 63), mid: (modulus.mid >> 1) | (modulus.hi << 63), hi: modulus.hi >> 1}
	return modulus, half, q01Lo, q01Hi
}

func q012Words(q0, q1, q2 uint64) (qLo, qMid, qHi, halfLo, halfMid, halfHi uint64) {
	modulus, half, _, _ := q012Modulus(q0, q1, q2)
	return modulus.lo, modulus.mid, modulus.hi, half.lo, half.mid, half.hi
}

func mod128By64(lo, hi, modulus uint64) uint64 {
	_, remainder := bits.Div64(0, hi%modulus, modulus)
	_, remainder = bits.Div64(remainder, lo, modulus)
	return remainder
}

func mod192By64(value uint192, modulus uint64) uint64 {
	_, remainder := bits.Div64(0, value.hi%modulus, modulus)
	_, remainder = bits.Div64(remainder, value.mid, modulus)
	_, remainder = bits.Div64(remainder, value.lo, modulus)
	return remainder
}

func crtQ012(r0, r1, r2, q0, q1, q2, q0InvQ1, q01InvQ2 uint64) uint192 {
	x01Lo, x01Hi := crtQ01(r0, r1, q0, q1, q0InvQ1)
	x01ModQ2 := mod128By64(x01Lo, x01Hi, q2)
	var delta uint64
	if r2 >= x01ModQ2 {
		delta = r2 - x01ModQ2
	} else {
		delta = q2 - (x01ModQ2 - r2)
	}
	prodHi, prodLo := bits.Mul64(delta, q01InvQ2)
	_, t := bits.Div64(prodHi, prodLo, q2)
	q01Hi, q01Lo := bits.Mul64(q0, q1)
	return add128To192(x01Lo, x01Hi, mul128By64(q01Lo, q01Hi, t))
}

func centeredQ012(value, modulus, half uint192) (uint192, bool) {
	if cmp192(value, half) > 0 {
		return sub192(modulus, value), true
	}
	return value, false
}

func roundedMagnitude192(value uint192, divisor uint64) uint192 {
	q2 := value.hi / divisor
	remainder := value.hi % divisor
	q1, remainder := bits.Div64(remainder, value.mid, divisor)
	q0, remainder := bits.Div64(remainder, value.lo, divisor)
	quotient := uint192{lo: q0, mid: q1, hi: q2}
	if remainder > divisor/2 {
		var carry uint64
		quotient.lo, carry = bits.Add64(quotient.lo, 1, 0)
		quotient.mid, carry = bits.Add64(quotient.mid, 0, carry)
		quotient.hi, _ = bits.Add64(quotient.hi, 0, carry)
	}
	return quotient
}

func signedResidue192(magnitude uint192, negative bool, modulus uint64) uint64 {
	remainder := mod192By64(magnitude, modulus)
	if negative && remainder != 0 {
		return modulus - remainder
	}
	return remainder
}
