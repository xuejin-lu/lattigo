package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// Evaluator is an explicit Fast CKKS execution boundary. It intentionally
// does not replace or modify ckks.Evaluator: callers select Fast execution by
// constructing this evaluator explicitly. Evaluators reuse internal scratch
// and are intended for one execution stream at a time; use separate instances
// for concurrent evaluation.
type Evaluator struct {
	Parameters     ckks.Parameters
	rescaleScratch fastRescaleScratch
	nttScratch     [3]ring.Poly
}

// NewEvaluator creates an explicit q0/q1-authoritative evaluator.
func NewEvaluator(params ckks.Parameters) *Evaluator {
	eval := &Evaluator{Parameters: params, rescaleScratch: newFastRescaleScratch(params.RingQ())}
	for i := range eval.nttScratch {
		eval.nttScratch[i] = ring.NewPoly(params.N(), 1)
	}
	return eval
}

// GetParameters returns the Fast evaluator's CKKS parameters.
func (eval *Evaluator) GetParameters() *ckks.Parameters { return &eval.Parameters }

// GetRLWEParameters satisfies rlwe.ParameterProvider without exposing a
// Standard evaluator or any full-Q/QP key-switching operation.
func (eval *Evaluator) GetRLWEParameters() *rlwe.Parameters { return &eval.Parameters.Parameters }

// Automorphism applies a Fast automorphism to a degree-one ciphertext. It
// operates only on q0 and q1 and does not use evaluation keys or relinearize.
// The input and output may alias. The ciphertexts must have the same level;
// Fast does not perform level reduction as part of this operation.
func (eval *Evaluator) Automorphism(ctIn, ctOut *rlwe.Ciphertext, galEl uint64) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if ctIn == nil || ctOut == nil {
		return errors.New("ctIn and ctOut cannot be nil")
	}
	if ctIn.MetaData == nil || ctOut.MetaData == nil {
		return errors.New("ctIn and ctOut metadata cannot be nil")
	}
	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		return errors.New("Fast automorphism requires degree-one ciphertexts")
	}
	if eval.Parameters.RingType() != ring.Standard {
		return fmt.Errorf("Fast automorphism requires the Standard ring, got %s", eval.Parameters.RingType())
	}
	if ctIn.Level() != ctOut.Level() || ctIn.Level() < 1 {
		return errors.New("Fast automorphism requires equal levels containing q0 and q1")
	}
	if ctIn.N() != eval.Parameters.N() || ctOut.N() != eval.Parameters.N() {
		return errors.New("ciphertext dimensions do not match Fast evaluator parameters")
	}
	if ctIn.IsNTT != ctOut.IsNTT || ctIn.IsMontgomery != ctOut.IsMontgomery {
		return errors.New("Fast automorphism requires equal input/output domains and Montgomery representations")
	}
	for name, ct := range map[string]*rlwe.Ciphertext{"ctIn": ctIn, "ctOut": ctOut} {
		for d := 0; d <= 1; d++ {
			if len(ct.Value[d].Coeffs) < 2 || len(ct.Value[d].Coeffs[0]) != eval.Parameters.N() || len(ct.Value[d].Coeffs[1]) != eval.Parameters.N() {
				return fmt.Errorf("%s has invalid q0/q1 storage", name)
			}
		}
	}

	ringQ := eval.Parameters.RingQ().AtLevel(ctIn.Level())
	// Keep the source intact until both authoritative components have been
	// transformed. This also makes ctIn == ctOut safe.
	tmp0 := ring.NewPoly(eval.Parameters.N(), 1)
	tmp1 := ring.NewPoly(eval.Parameters.N(), 1)
	if err := FastAutomorphism(ringQ, ctIn.Value[0], tmp0, galEl, ctIn.IsNTT); err != nil {
		return fmt.Errorf("FastAutomorphism(c0): %w", err)
	}
	if err := FastAutomorphism(ringQ, ctIn.Value[1], tmp1, galEl, ctIn.IsNTT); err != nil {
		return fmt.Errorf("FastAutomorphism(c1): %w", err)
	}
	copyQ01(tmp0, ctOut.Value[0])
	copyQ01(tmp1, ctOut.Value[1])
	*ctOut.MetaData = *ctIn.MetaData
	return nil
}

// Rotate applies the Fast automorphism corresponding to a CKKS slot rotation.
func (eval *Evaluator) Rotate(ctIn, ctOut *rlwe.Ciphertext, k int) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	return eval.Automorphism(ctIn, ctOut, eval.Parameters.GaloisElementForRotation(k))
}

// RotateNew returns a newly allocated ciphertext containing the Fast rotation
// of ctIn.
func (eval *Evaluator) RotateNew(ctIn *rlwe.Ciphertext, k int) (*rlwe.Ciphertext, error) {
	if eval == nil || ctIn == nil {
		return nil, errors.New("evaluator and ciphertext cannot be nil")
	}
	ctOut := ckks.NewCiphertext(eval.Parameters, 1, ctIn.Level())
	ctOut.IsNTT = ctIn.IsNTT
	ctOut.IsMontgomery = ctIn.IsMontgomery
	return ctOut, eval.Rotate(ctIn, ctOut, k)
}

// Conjugate applies the CKKS complex-conjugation automorphism without an
// evaluation key. Fast supports this only for the Standard ring, through the
// same q0/q1 automorphism path as Rotate.
func (eval *Evaluator) Conjugate(ctIn, ctOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if eval.Parameters.RingType() != ring.Standard {
		return fmt.Errorf("Fast Conjugate requires the Standard ring, got %s", eval.Parameters.RingType())
	}
	return eval.Automorphism(ctIn, ctOut, eval.Parameters.GaloisElementForComplexConjugation())
}

// ConjugateNew returns the Fast conjugate of ctIn.
func (eval *Evaluator) ConjugateNew(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil || ctIn == nil {
		return nil, errors.New("evaluator and ciphertext cannot be nil")
	}
	ctOut := ckks.NewCiphertext(eval.Parameters, 1, ctIn.Level())
	ctOut.IsNTT = ctIn.IsNTT
	ctOut.IsMontgomery = ctIn.IsMontgomery
	return ctOut, eval.Conjugate(ctIn, ctOut)
}
