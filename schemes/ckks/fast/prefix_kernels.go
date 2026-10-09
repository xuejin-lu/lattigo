package fast

import (
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

// validatePrefixRows validates an explicitly requested Q-prefix operation.
// The row count is never inferred from allocated backing: callers must choose
// the authoritative width they intend to process.
func validatePrefixRows(ringQ *ring.Ring, level, rows int, polys ...ring.Poly) error {
	return fastcore.ValidatePrefixRows(ringQ, level, rows, polys...)
}

// validatePrefixBacking checks physical scratch/operand row availability when
// the buffer itself does not represent a ciphertext logical Level (for
// example, evaluator scratch stores only the capped four-row prefix).
func validatePrefixBacking(ringQ *ring.Ring, rows int, polys ...ring.Poly) error {
	return fastcore.ValidatePrefixBacking(ringQ, rows, polys...)
}

func addSubPrefixRows(ringQ *ring.Ring, level, rows int, a, b, out ring.Poly, sub bool) error {
	return fastcore.AddSubPolynomialRows(ringQ, level, rows, sub, fastcore.AddSubPolynomialPair{Op0: a, Op1: b, Output: out})
}

func copyPrefixRows(ringQ *ring.Ring, level, rows int, src, dst ring.Poly) error {
	return fastcore.CopyPrefixRows(ringQ, level, rows, src, dst)
}

func copyPrefixRowsUnchecked(rows int, src, dst ring.Poly) {
	fastcore.CopyPrefixRowsUnchecked(rows, src, dst)
}

func negatePrefixRows(ringQ *ring.Ring, level, rows int, src, dst ring.Poly) error {
	return fastcore.NegatePrefixRows(ringQ, level, rows, src, dst)
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
