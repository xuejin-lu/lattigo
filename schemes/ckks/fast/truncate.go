package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// FastTruncateDegree2To1 discards the degree-two component of a Fast
// ciphertext. Under Fast secret semantics (s = 0), c2 does not contribute to
// the represented value. Only the authoritative q0 and q1 limbs of c0 and c1
// are copied; dormant limbs are deliberately neither read nor written.
func FastTruncateDegree2To1(ctIn, ctOut *rlwe.Ciphertext) error {
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
	if ctIn.Level() < 1 {
		return errors.New("FastTruncateDegree2To1 requires at least q0 and q1")
	}
	for i := 0; i < 2; i++ {
		if len(ctIn.Value[i].Coeffs) < 2 || len(ctOut.Value[i].Coeffs) < 2 ||
			len(ctIn.Value[i].Coeffs[0]) != ctIn.N() || len(ctIn.Value[i].Coeffs[1]) != ctIn.N() ||
			len(ctOut.Value[i].Coeffs[0]) != ctOut.N() || len(ctOut.Value[i].Coeffs[1]) != ctOut.N() {
			return fmt.Errorf("ciphertext component %d has invalid q0/q1 storage", i)
		}
	}

	// Resize only drops the c2 component. It does not inspect any polynomial
	// coefficients, which makes this safe when ctIn and ctOut alias.
	ctOut.Resize(1, ctIn.Level())
	if ctOut != ctIn {
		copyQ01(ctIn.Value[0], ctOut.Value[0])
		copyQ01(ctIn.Value[1], ctOut.Value[1])
	}
	*ctOut.MetaData = *ctIn.MetaData
	return nil
}
