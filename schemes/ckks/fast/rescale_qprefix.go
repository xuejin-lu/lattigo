package fast

import (
	"errors"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// rescaleNQPrefix is a thin facade over the same import-neutral transactional
// core used by the ordinary public CKKS evaluator.
func (eval *Evaluator) rescaleNQPrefix(op0 *rlwe.Ciphertext, nbRescales, sourceRows int, opOut *rlwe.Ciphertext) error {
	if eval == nil || eval.rescaleCore == nil {
		return errors.New("Fast Rescale core is not initialized")
	}
	return eval.rescaleCore.ApplyRows(eval.Parameters.RingQ(), op0, opOut, nbRescales, sourceRows)
}

// The following exact coefficient/reconstruction helpers remain for the
// existing capacity observer, scalar guard, and level-transition diagnostics.
// Production Rescale itself is implemented only by fastcore.RescaleWorkspace.
func prefixToCoefficientRows(ringQ *ring.Ring, src ring.Poly, rows int, isNTT, isMontgomery bool, dst ring.Poly) error {
	if err := validatePrefixBacking(ringQ, rows, src, dst); err != nil {
		return err
	}
	for row := 0; row < rows; row++ {
		if isNTT {
			ringQ.SubRings[row].INTT(src.Coeffs[row], dst.Coeffs[row])
		} else {
			copy(dst.Coeffs[row], src.Coeffs[row])
		}
		if isMontgomery {
			ringQ.SubRings[row].IMForm(dst.Coeffs[row], dst.Coeffs[row])
		}
	}
	return nil
}

func reconstructQPrefix(rows int, residues [MaxQPrefixWidth]uint64, scratch *fastRescaleScratch) uint192 {
	switch rows {
	case 1:
		return uint192{lo: residues[0]}
	case 2:
		lo, hi := crtQ01(residues[0], residues[1], scratch.q[0], scratch.q[1], scratch.inverse[1])
		return uint192{lo: lo, mid: hi}
	case 3:
		return crtQ012Prepared(residues[0], residues[1], residues[2], scratch.q[0], scratch.q[1], scratch.q[2], scratch.inverse[1], scratch.inverse[2], scratch.modulus[1])
	case 4:
		return crtQ0123Prepared(residues[0], residues[1], residues[2], residues[3], scratch.q[0], scratch.q[1], scratch.q[2], scratch.q[3], scratch.inverse[1], scratch.inverse[2], scratch.inverse[3], scratch.modulus[1], scratch.modulus[2])
	default:
		return uint192{}
	}
}

func centeredQPrefix(value, modulus, half uint192) (magnitude uint192, negative bool) {
	if cmp192(value, half) > 0 {
		return sub192(modulus, value), true
	}
	return value, false
}

func uint192ToBigInt(value uint192) *big.Int {
	out := new(big.Int).SetUint64(value.hi)
	out.Lsh(out, 64)
	out.Add(out, new(big.Int).SetUint64(value.mid))
	out.Lsh(out, 64)
	out.Add(out, new(big.Int).SetUint64(value.lo))
	return out
}
