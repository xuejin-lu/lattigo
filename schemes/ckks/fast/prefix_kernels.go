package fast

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// validatePrefixRows validates an explicitly requested Q-prefix operation.
// The row count is never inferred from allocated backing: callers must choose
// the authoritative width they intend to process.
func validatePrefixRows(ringQ *ring.Ring, level, rows int, polys ...ring.Poly) error {
	if ringQ == nil {
		return fmt.Errorf("Q-prefix ring cannot be nil")
	}
	if level < 0 || level > ringQ.Level() {
		return fmt.Errorf("Q-prefix level %d is outside ring level [0,%d]", level, ringQ.Level())
	}
	width, err := QPrefixWidth(level)
	if err != nil {
		return err
	}
	if rows < 1 || rows > level+1 || rows > width {
		return fmt.Errorf("Q-prefix row count %d must be in [1,%d] at level %d", rows, min(level+1, width), level)
	}
	if rows > len(ringQ.SubRings) {
		return fmt.Errorf("Q-prefix row count %d exceeds available q subrings %d", rows, len(ringQ.SubRings))
	}
	for row := 0; row < rows; row++ {
		if ringQ.SubRings[row] == nil {
			return fmt.Errorf("q%d subring is unavailable", row)
		}
	}
	for i, poly := range polys {
		if len(poly.Coeffs) == 0 || poly.Level() < level {
			return fmt.Errorf("polynomial %d logical level is below requested level %d", i, level)
		}
	}
	return validatePrefixBacking(ringQ, rows, polys...)
}

// validatePrefixBacking checks physical scratch/operand row availability when
// the buffer itself does not represent a ciphertext logical Level (for
// example, evaluator scratch stores only the capped four-row prefix).
func validatePrefixBacking(ringQ *ring.Ring, rows int, polys ...ring.Poly) error {
	if ringQ == nil {
		return fmt.Errorf("Q-prefix ring cannot be nil")
	}
	if rows < 1 || rows > MaxQPrefixWidth || rows > len(ringQ.SubRings) {
		return fmt.Errorf("Q-prefix backing row count %d is unavailable", rows)
	}
	for row := 0; row < rows; row++ {
		if ringQ.SubRings[row] == nil {
			return fmt.Errorf("q%d subring is unavailable", row)
		}
	}
	for i, poly := range polys {
		if len(poly.Coeffs) == 0 {
			return fmt.Errorf("polynomial %d has no coefficient rows", i)
		}
		if poly.N() != ringQ.N() {
			return fmt.Errorf("polynomial %d dimension %d does not match ring dimension %d", i, poly.N(), ringQ.N())
		}
		if len(poly.Coeffs) < rows {
			return fmt.Errorf("polynomial %d has %d rows, requested %d", i, len(poly.Coeffs), rows)
		}
		for row := 0; row < rows; row++ {
			if len(poly.Coeffs[row]) != ringQ.N() {
				return fmt.Errorf("polynomial %d q%d backing length %d does not match N=%d", i, row, len(poly.Coeffs[row]), ringQ.N())
			}
		}
	}
	return nil
}

func addSubPrefixRows(ringQ *ring.Ring, level, rows int, a, b, out ring.Poly, sub bool) error {
	if err := validatePrefixRows(ringQ, level, rows, a, b, out); err != nil {
		return err
	}
	for row := 0; row < rows; row++ {
		if sub {
			ringQ.SubRings[row].Sub(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		} else {
			ringQ.SubRings[row].Add(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		}
	}
	return nil
}

func copyPrefixRows(ringQ *ring.Ring, level, rows int, src, dst ring.Poly) error {
	if err := validatePrefixRows(ringQ, level, rows, src, dst); err != nil {
		return err
	}
	copyPrefixRowsUnchecked(rows, src, dst)
	return nil
}

func copyPrefixRowsUnchecked(rows int, src, dst ring.Poly) {
	for row := 0; row < rows; row++ {
		copy(dst.Coeffs[row], src.Coeffs[row])
	}
}

func negatePrefixRows(ringQ *ring.Ring, level, rows int, src, dst ring.Poly) error {
	if err := validatePrefixRows(ringQ, level, rows, src, dst); err != nil {
		return err
	}
	for row := 0; row < rows; row++ {
		ringQ.SubRings[row].Neg(src.Coeffs[row], dst.Coeffs[row])
	}
	return nil
}

func pointMulPrefixRows(ringQ *ring.Ring, level, rows int, a, b, out ring.Poly, montgomery, add bool) error {
	if err := validatePrefixRows(ringQ, level, rows, a, b); err != nil {
		return err
	}
	// The output may be a compact evaluator scratch polynomial whose header
	// level is only its physical prefix width, not the ciphertext operation's
	// logical level. The explicit level/rows pair governs the operation; still
	// require every output row to have N-sized backing before writing.
	if err := validatePrefixBacking(ringQ, rows, out); err != nil {
		return err
	}
	for row := 0; row < rows; row++ {
		subring := ringQ.SubRings[row]
		switch {
		case montgomery && add:
			subring.MulCoeffsMontgomeryThenAdd(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		case montgomery:
			subring.MulCoeffsMontgomery(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		case add:
			subring.MulCoeffsBarrettThenAdd(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		default:
			subring.MulCoeffsBarrett(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		}
	}
	return nil
}
