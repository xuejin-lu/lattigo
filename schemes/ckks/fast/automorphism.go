package fast

import (
	"errors"
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastAutomorphism applies a ring automorphism using the legacy authoritative
// row count. Higher rows are deliberately neither read nor written. The output
// may alias the input.
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
	if ringQ.Level() < 1 || polIn.Level() < 1 || polOut.Level() < 1 {
		return errors.New("Fast automorphism requires q0 and q1")
	}
	level := min(ringQ.Level(), minPolyLevel(polIn, polOut))
	rows := maintainedLimbCountForRingAtLevel(ringQ, level)
	return FastAutomorphismRows(ringQ, polIn, polOut, galEl, isNTT, rows)
}

// FastAutomorphismRows applies a ring automorphism to exactly rows explicitly
// requested Q-prefix rows. Higher rows are neither read nor written. The output
// may alias the input.
func FastAutomorphismRows(ringQ *ring.Ring, polIn, polOut ring.Poly, galEl uint64, isNTT bool, rows int) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Type() != ring.Standard {
		return fmt.Errorf("Fast automorphism requires the Standard ring, got %s", ringQ.Type())
	}
	level := min(ringQ.Level(), minPolyLevel(polIn, polOut))
	if err := validatePrefixRows(ringQ, level, rows, polIn, polOut); err != nil {
		return err
	}
	if isNTT {
		index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galEl)
		if err != nil {
			return fmt.Errorf("compute NTT automorphism index: %w", err)
		}
		for limb := 0; limb < rows; limb++ {
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
	for limb := 0; limb < rows; limb++ {
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

// fastAutomorphism is the evaluator-owned legacy-width hot path. Its index
// cache and scratch buffers are deliberately scoped to one Evaluator.
func (eval *Evaluator) fastAutomorphism(ringQ *ring.Ring, polIn, polOut ring.Poly, galEl uint64, isNTT bool) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Level() < 1 || polIn.Level() < 1 || polOut.Level() < 1 {
		return errors.New("Fast automorphism requires q0 and q1")
	}
	level := min(ringQ.Level(), minPolyLevel(polIn, polOut))
	rows := maintainedLimbCountForRingAtLevel(ringQ, level)
	return eval.fastAutomorphismRows(ringQ, polIn, polOut, galEl, isNTT, rows)
}

func (eval *Evaluator) fastAutomorphismRows(ringQ *ring.Ring, polIn, polOut ring.Poly, galEl uint64, isNTT bool, rows int) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Type() != ring.Standard {
		return fmt.Errorf("Fast automorphism requires the Standard ring, got %s", ringQ.Type())
	}
	level := min(ringQ.Level(), minPolyLevel(polIn, polOut))
	if err := validatePrefixRows(ringQ, level, rows, polIn, polOut); err != nil {
		return err
	}
	if len(eval.automorphismScratch) < rows {
		return fmt.Errorf("automorphism scratch has %d rows, requested %d", len(eval.automorphismScratch), rows)
	}
	for row := 0; row < rows; row++ {
		if len(eval.automorphismScratch[row]) != ringQ.N() {
			return fmt.Errorf("automorphism scratch q%d has invalid dimension", row)
		}
	}
	if !isNTT {
		return FastAutomorphismRows(ringQ, polIn, polOut, galEl, false, rows)
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
	for limb := 0; limb < rows; limb++ {
		tmp := eval.automorphismScratch[limb]
		for j, src := range index {
			tmp[j] = polIn.Coeffs[limb][src]
		}
		copy(polOut.Coeffs[limb], tmp)
	}
	return nil
}
