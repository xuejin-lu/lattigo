package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastPartialNTT applies the forward NTT only to q0 and q1 of p1 and writes
// those two limbs to p2. All limbs with index >= 2 in p2 are left untouched.
//
// The operation is deliberately independent of rlwe metadata: ring.Poly does
// not carry domain state, so the caller must track the Fast mixed-domain state
// separately. Inputs and outputs are the same representations expected by
// ring.SubRing.NTT. In particular, this function does not perform a
// Montgomery conversion or update an IsNTT flag.
func FastPartialNTT(ringQ *ring.Ring, p1, p2 ring.Poly) error {
	if err := validatePartialTransform(ringQ, p1, p2); err != nil {
		return err
	}

	for limb := 0; limb < maintainedLimbCountForRingAtLevel(ringQ, minPolyLevel(p1, p2)); limb++ {
		ringQ.SubRings[limb].NTT(p1.Coeffs[limb], p2.Coeffs[limb])
	}
	return nil
}

// FastPartialINTT applies the inverse NTT only to q0 and q1 of p1 and writes
// those two limbs to p2. All limbs with index >= 2 in p2 are left untouched.
// The output of the underlying SubRing.INTT is canonical coefficient-domain
// residue data and is therefore suitable for the Phase 1A CRT precondition
// once the caller has also ensured that the values are not Montgomery data.
func FastPartialINTT(ringQ *ring.Ring, p1, p2 ring.Poly) error {
	if err := validatePartialTransform(ringQ, p1, p2); err != nil {
		return err
	}

	for limb := 0; limb < maintainedLimbCountForRingAtLevel(ringQ, minPolyLevel(p1, p2)); limb++ {
		ringQ.SubRings[limb].INTT(p1.Coeffs[limb], p2.Coeffs[limb])
	}
	return nil
}

func validatePartialTransform(ringQ *ring.Ring, p1, p2 ring.Poly) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Level() < 1 {
		return errors.New("FastPartialNTT/INTT requires q0 and q1")
	}
	if p1.Level() < 1 || p2.Level() < 1 {
		return errors.New("input and output require q0 and q1")
	}
	if p1.N() != ringQ.N() || p2.N() != ringQ.N() {
		return fmt.Errorf("input and output dimensions must match ringQ.N()=%d", ringQ.N())
	}
	maintained := maintainedLimbCountForRingAtLevel(ringQ, minPolyLevel(p1, p2))
	if len(p1.Coeffs) < maintained || len(p2.Coeffs) < maintained {
		return fmt.Errorf("maintained coefficient slices are insufficient")
	}
	for limb := 0; limb < maintained; limb++ {
		if len(p1.Coeffs[limb]) != ringQ.N() || len(p2.Coeffs[limb]) != ringQ.N() {
			return fmt.Errorf("maintained coefficient slices must have length %d", ringQ.N())
		}
	}
	return nil
}

func minPolyLevel(a, b ring.Poly) int {
	level := a.Level()
	if b.Level() < level {
		level = b.Level()
	}
	return level
}
