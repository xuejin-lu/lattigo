package fast

import (
	"errors"
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastAutomorphism applies a ring automorphism to the authoritative q0 and q1
// limbs of a polynomial. Limbs q2 and above are deliberately neither read nor
// written. The output may alias the input.
//
// Fast currently supports the Standard ring only. In particular, this helper
// does not implement the different folded-index semantics of the
// conjugate-invariant ring.
func FastAutomorphism(ringQ *ring.Ring, polIn, polOut ring.Poly, galEl uint64, isNTT bool) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Type() != ring.Standard {
		return fmt.Errorf("Fast automorphism requires the Standard ring, got %s", ringQ.Type())
	}
	if ringQ.Level() < 1 {
		return errors.New("Fast automorphism requires q0 and q1")
	}
	if len(polIn.Coeffs) < 2 || len(polOut.Coeffs) < 2 {
		return errors.New("Fast automorphism requires q0 and q1 polynomial limbs")
	}
	if len(polIn.Coeffs[0]) != ringQ.N() || len(polIn.Coeffs[1]) != ringQ.N() ||
		len(polOut.Coeffs[0]) != ringQ.N() || len(polOut.Coeffs[1]) != ringQ.N() {
		return errors.New("Fast automorphism polynomial dimensions do not match ringQ")
	}

	if isNTT {
		index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galEl)
		if err != nil {
			return fmt.Errorf("compute NTT automorphism index: %w", err)
		}
		for limb := 0; limb < 2; limb++ {
			// The temporary is intentional: ring automorphism primitives are
			// non-in-place, while Fast explicitly permits input/output aliasing.
			tmp := make([]uint64, ringQ.N())
			for j, src := range index {
				tmp[j] = polIn.Coeffs[limb][src]
			}
			copy(polOut.Coeffs[limb], tmp)
		}
		return nil
	}

	mask := uint64(ringQ.N() - 1)
	logN := uint(bits.Len64(mask))
	for limb := 0; limb < 2; limb++ {
		modulus := ringQ.SubRings[limb].Modulus
		tmp := make([]uint64, ringQ.N())
		for i, value := range polIn.Coeffs[limb] {
			raw := uint64(i) * galEl
			index := raw & mask
			if (raw>>logN)&1 == 0 {
				tmp[index] = value
			} else {
				// Keep the same representation as ring.Ring.Automorphism,
				// including modulus for an input zero.
				tmp[index] = modulus - value
			}
		}
		copy(polOut.Coeffs[limb], tmp)
	}
	return nil
}

// fastAutomorphism is the evaluator-owned hot path. Its index cache and
// q0/q1 scratch buffers are deliberately scoped to one Evaluator, which is
// already documented as a single execution stream.
func (eval *Evaluator) fastAutomorphism(ringQ *ring.Ring, polIn, polOut ring.Poly, galEl uint64, isNTT bool) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Type() != ring.Standard {
		return fmt.Errorf("Fast automorphism requires the Standard ring, got %s", ringQ.Type())
	}
	if ringQ.Level() < 1 {
		return errors.New("Fast automorphism requires q0 and q1")
	}
	if len(polIn.Coeffs) < 2 || len(polOut.Coeffs) < 2 {
		return errors.New("Fast automorphism requires q0 and q1 polynomial limbs")
	}
	if len(polIn.Coeffs[0]) != ringQ.N() || len(polIn.Coeffs[1]) != ringQ.N() ||
		len(polOut.Coeffs[0]) != ringQ.N() || len(polOut.Coeffs[1]) != ringQ.N() {
		return errors.New("Fast automorphism polynomial dimensions do not match ringQ")
	}
	if !isNTT {
		return FastAutomorphism(ringQ, polIn, polOut, galEl, false)
	}

	index, ok := eval.automorphismIndexCache[galEl]
	if !ok {
		var err error
		index, err = ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galEl)
		if err != nil {
			return fmt.Errorf("compute NTT automorphism index: %w", err)
		}
		eval.automorphismIndexCache[galEl] = index
	}
	for limb := 0; limb < 2; limb++ {
		tmp := eval.automorphismScratch[limb]
		for j, src := range index {
			tmp[j] = polIn.Coeffs[limb][src]
		}
		copy(polOut.Coeffs[limb], tmp)
	}
	return nil
}
