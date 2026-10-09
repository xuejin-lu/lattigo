package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
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
	workspace := fastcore.NewAutomorphismWorkspace(ringQ.N(), rows)
	return workspace.ApplyRows(ringQ, galEl, isNTT, rows, fastcore.PolynomialPair{Input: polIn, Output: polOut})
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
	if eval.automorphismCore == nil {
		return errors.New("Fast automorphism core is not initialized")
	}
	return eval.automorphismCore.ApplyRows(ringQ, galEl, isNTT, rows,
		fastcore.PolynomialPair{Input: polIn, Output: polOut})
}
