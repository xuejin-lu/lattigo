package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastTruncateDegree2To1 discards the degree-two component of a Fast
// ciphertext. Under Fast secret semantics (s = 0), c2 does not contribute to
// the represented value. Only the legacy authoritative rows of c0 and c1 are
// copied; higher rows are deliberately neither read nor written.
func FastTruncateDegree2To1(ctIn, ctOut *rlwe.Ciphertext) error {
	return fastTruncateDegree2To1(nil, ctIn, ctOut)
}

// fastTruncateDegree2To1 is the evaluator-internal legacy-width variant. The
// public helper above intentionally retains its legacy q0/q1-only semantics.
func fastTruncateDegree2To1(ringQ *ring.Ring, ctIn, ctOut *rlwe.Ciphertext) error {
	maintained := 2
	if ringQ != nil && ctIn != nil {
		maintained = maintainedLimbCountForRingAtLevel(ringQ, ctIn.Level())
	}
	return fastTruncateDegree2To1Rows(ringQ, ctIn, ctOut, maintained)
}

// fastTruncateDegree2To1Rows copies exactly rows explicitly requested
// Q-prefix rows from c0/c1 and discards c2. It supports Level 0 when rows is 1.
func fastTruncateDegree2To1Rows(ringQ *ring.Ring, ctIn, ctOut *rlwe.Ciphertext, rows int) error {
	if ctIn == nil || ctOut == nil {
		return errors.New("input and output ciphertexts cannot be nil")
	}
	if ctIn.MetaData == nil || ctOut.MetaData == nil {
		return errors.New("input and output ciphertext metadata cannot be nil")
	}
	if ctIn.Degree() != 2 {
		return fmt.Errorf("FastTruncateDegree2To1 requires degree-two input, got degree %d", ctIn.Degree())
	}
	if ctOut.Degree() < 1 {
		return fmt.Errorf("FastTruncateDegree2To1 requires output storage for degree one, got degree %d", ctOut.Degree())
	}
	if ctIn.N() != ctOut.N() {
		return errors.New("input and output ciphertext dimensions do not match")
	}
	if ctIn.Level() != ctOut.Level() {
		return fmt.Errorf("input and output ciphertext levels do not match: %d != %d", ctIn.Level(), ctOut.Level())
	}
	if ringQ == nil && ctIn.Level() < 1 {
		return errors.New("FastTruncateDegree2To1 requires at least q0 and q1")
	}
	if ringQ == nil {
		if rows != 2 {
			return errors.New("legacy FastTruncateDegree2To1 requires exactly q0/q1 rows")
		}
	} else if err := validatePrefixRows(ringQ, ctIn.Level(), rows, ctIn.Value[0], ctOut.Value[0], ctIn.Value[1], ctOut.Value[1]); err != nil {
		return err
	}
	for i := 0; i < 2; i++ {
		if ringQ == nil {
			if len(ctIn.Value[i].Coeffs) < rows || len(ctOut.Value[i].Coeffs) < rows {
				return fmt.Errorf("ciphertext component %d has insufficient maintained-limb storage", i)
			}
		}
		for limb := 0; ringQ == nil && limb < rows; limb++ {
			if len(ctIn.Value[i].Coeffs[limb]) != ctIn.N() || len(ctOut.Value[i].Coeffs[limb]) != ctOut.N() {
				return fmt.Errorf("ciphertext component %d has invalid maintained limb %d storage", i, limb)
			}
		}
	}

	// Resize only drops the c2 component. It does not inspect any polynomial
	// coefficients, which makes this safe when ctIn and ctOut alias.
	Resize(ctOut, 1, ctIn.Level(), ctIn.N())
	if ctOut != ctIn {
		copyPrefixRowsUnchecked(rows, ctIn.Value[0], ctOut.Value[0])
		copyPrefixRowsUnchecked(rows, ctIn.Value[1], ctOut.Value[1])
	}
	*ctOut.MetaData = *ctIn.MetaData
	return nil
}
