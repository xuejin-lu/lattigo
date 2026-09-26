package fast

import (
	"errors"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastPartialNTT applies the forward NTT using the legacy authoritative row
// count. All higher rows in p2 are left untouched.
//
// The operation is deliberately independent of rlwe metadata: ring.Poly does
// not carry domain state, so the caller must track the Fast mixed-domain state
// separately. Inputs and outputs are the same representations expected by
// ring.SubRing.NTT. In particular, this function does not perform a
// Montgomery conversion or update an IsNTT flag.
func FastPartialNTT(ringQ *ring.Ring, p1, p2 ring.Poly) error {
	if ringQ == nil || ringQ.Level() < 1 || p1.Level() < 1 || p2.Level() < 1 {
		return errors.New("FastPartialNTT requires q0 and q1")
	}
	level := min(ringQ.Level(), minPolyLevel(p1, p2))
	rows := maintainedLimbCountForRingAtLevel(ringQ, level)
	return FastPartialNTTRows(ringQ, p1, p2, rows)
}

// FastPartialNTTRows applies NTT to exactly rows explicitly requested prefix
// rows. Rows outside that count are not read or written.
func FastPartialNTTRows(ringQ *ring.Ring, p1, p2 ring.Poly, rows int) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	level := min(ringQ.Level(), minPolyLevel(p1, p2))
	if err := validatePrefixRows(ringQ, level, rows, p1, p2); err != nil {
		return err
	}
	for limb := 0; limb < rows; limb++ {
		ringQ.SubRings[limb].NTT(p1.Coeffs[limb], p2.Coeffs[limb])
	}
	return nil
}

// FastPartialINTT applies the inverse NTT using the legacy authoritative row
// count. All higher rows in p2 are left untouched.
// The output of the underlying SubRing.INTT is canonical coefficient-domain
// residue data and is therefore suitable for the Phase 1A CRT precondition
// once the caller has also ensured that the values are not Montgomery data.
func FastPartialINTT(ringQ *ring.Ring, p1, p2 ring.Poly) error {
	if ringQ == nil || ringQ.Level() < 1 || p1.Level() < 1 || p2.Level() < 1 {
		return errors.New("FastPartialINTT requires q0 and q1")
	}
	level := min(ringQ.Level(), minPolyLevel(p1, p2))
	rows := maintainedLimbCountForRingAtLevel(ringQ, level)
	return FastPartialINTTRows(ringQ, p1, p2, rows)
}

// FastPartialINTTRows applies INTT to exactly rows explicitly requested
// prefix rows. Rows outside that count are not read or written.
func FastPartialINTTRows(ringQ *ring.Ring, p1, p2 ring.Poly, rows int) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	level := min(ringQ.Level(), minPolyLevel(p1, p2))
	if err := validatePrefixRows(ringQ, level, rows, p1, p2); err != nil {
		return err
	}
	for limb := 0; limb < rows; limb++ {
		ringQ.SubRings[limb].INTT(p1.Coeffs[limb], p2.Coeffs[limb])
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
