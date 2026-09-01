package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// LinearTransform evaluates a single-level diagonal linear transformation
// using only the authoritative q0 and q1 limbs. Each diagonal is evaluated as
// FastAutomorphism followed by q0/q1 plaintext multiplication and accumulation.
// No evaluation key, QP basis, or level transition is involved.
//
// The transformation must have the same level, domain, Montgomery
// representation, and ring degree as the input ciphertext. The output may
// alias the input; dormant limbs in the output are left untouched.
func (eval *Evaluator) LinearTransform(ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation, ctOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if ctIn == nil || ctOut == nil {
		return errors.New("ctIn and ctOut cannot be nil")
	}
	if ctIn.MetaData == nil || ctOut.MetaData == nil || matrix.MetaData == nil {
		return errors.New("ciphertext and linear transformation metadata cannot be nil")
	}
	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		return errors.New("Fast LinearTransform requires degree-one ciphertexts")
	}
	if eval.Parameters.RingType() != ring.Standard {
		return fmt.Errorf("Fast LinearTransform requires the Standard ring, got %s", eval.Parameters.RingType())
	}
	if ctIn.N() != eval.Parameters.N() || ctOut.N() != eval.Parameters.N() {
		return errors.New("ciphertext dimensions do not match Fast evaluator parameters")
	}
	if ctIn.Level() != ctOut.Level() || ctIn.Level() < 1 {
		return errors.New("Fast LinearTransform requires equal levels containing q0 and q1")
	}
	if matrix.LevelQ != ctIn.Level() {
		return errors.New("Fast LinearTransform requires matrix.LevelQ to equal the ciphertext level")
	}
	if !ctIn.IsNTT || !matrix.IsNTT || !ctIn.IsMontgomery || !matrix.IsMontgomery {
		return errors.New("Fast LinearTransform requires NTT and Montgomery representation")
	}
	if ctIn.IsBatched != matrix.IsBatched {
		return errors.New("Fast LinearTransform requires matching batching metadata")
	}
	if len(ctIn.Value) < 2 || len(ctOut.Value) < 2 || len(ctIn.Value[0].Coeffs) < 2 || len(ctIn.Value[1].Coeffs) < 2 || len(ctOut.Value[0].Coeffs) < 2 || len(ctOut.Value[1].Coeffs) < 2 {
		return errors.New("Fast LinearTransform requires q0/q1 ciphertext storage")
	}

	ringQ := eval.Parameters.RingQ().AtLevel(ctIn.Level())
	acc0 := ring.NewPoly(eval.Parameters.N(), 1)
	acc1 := ring.NewPoly(eval.Parameters.N(), 1)
	rot0 := ring.NewPoly(eval.Parameters.N(), 1)
	rot1 := ring.NewPoly(eval.Parameters.N(), 1)
	term0 := ring.NewPoly(eval.Parameters.N(), 1)
	term1 := ring.NewPoly(eval.Parameters.N(), 1)

	first := true
	slots := 1 << matrix.LogDimensions.Cols
	for diagonal, plaintext := range matrix.Vec {
		diagonal &= slots - 1
		if diagonal >= slots {
			return fmt.Errorf("diagonal %d is outside the matrix slot range", diagonal)
		}
		if err := validateFastDiagonal(ringQ, plaintext.Q); err != nil {
			return fmt.Errorf("diagonal %d: %w", diagonal, err)
		}

		galEl := eval.Parameters.GaloisElement(diagonal)
		if err := FastAutomorphism(ringQ, ctIn.Value[0], rot0, galEl, ctIn.IsNTT); err != nil {
			return fmt.Errorf("diagonal %d automorphism(c0): %w", diagonal, err)
		}
		if err := FastAutomorphism(ringQ, ctIn.Value[1], rot1, galEl, ctIn.IsNTT); err != nil {
			return fmt.Errorf("diagonal %d automorphism(c1): %w", diagonal, err)
		}
		if err := fastPlaintextMul(ringQ, plaintext.Q, rot0, term0); err != nil {
			return fmt.Errorf("diagonal %d plaintext multiplication(c0): %w", diagonal, err)
		}
		if err := fastPlaintextMul(ringQ, plaintext.Q, rot1, term1); err != nil {
			return fmt.Errorf("diagonal %d plaintext multiplication(c1): %w", diagonal, err)
		}

		if first {
			copyQ01(term0, acc0)
			copyQ01(term1, acc1)
			first = false
		} else {
			for limb := 0; limb < 2; limb++ {
				ringQ.SubRings[limb].Add(term0.Coeffs[limb], acc0.Coeffs[limb], acc0.Coeffs[limb])
				ringQ.SubRings[limb].Add(term1.Coeffs[limb], acc1.Coeffs[limb], acc1.Coeffs[limb])
			}
		}
	}

	if first {
		return errors.New("Fast LinearTransform requires at least one diagonal")
	}

	*ctOut.MetaData = *ctIn.MetaData
	ctOut.Scale = ctIn.Scale.Mul(matrix.Scale)
	copyQ01(acc0, ctOut.Value[0])
	copyQ01(acc1, ctOut.Value[1])
	return nil
}

// fastPlaintextMul multiplies a q0/q1 plaintext diagonal by a q0/q1
// ciphertext polynomial. Both operands are in NTT and Montgomery form. Only
// SubRing operations are used.
func fastPlaintextMul(ringQ *ring.Ring, plaintext, ciphertext, output ring.Poly) error {
	if err := validateFastDiagonal(ringQ, plaintext); err != nil {
		return err
	}
	if ciphertext.N() != ringQ.N() || output.N() != ringQ.N() || len(ciphertext.Coeffs) < 2 || len(output.Coeffs) < 2 {
		return errors.New("ciphertext/output dimensions do not match ringQ")
	}
	for limb := 0; limb < 2; limb++ {
		ringQ.SubRings[limb].MulCoeffsMontgomery(plaintext.Coeffs[limb], ciphertext.Coeffs[limb], output.Coeffs[limb])
	}
	return nil
}

func validateFastDiagonal(ringQ *ring.Ring, plaintext ring.Poly) error {
	if plaintext.N() != ringQ.N() || plaintext.Level() < 1 {
		return errors.New("diagonal dimensions or level are insufficient")
	}
	if len(plaintext.Coeffs) < 2 || len(plaintext.Coeffs[0]) != ringQ.N() || len(plaintext.Coeffs[1]) != ringQ.N() {
		return errors.New("diagonal q0/q1 storage is invalid")
	}
	return nil
}

// FastLinearTransform is the functional form of Evaluator.LinearTransform.
func FastLinearTransform(params ckks.Parameters, ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation, ctOut *rlwe.Ciphertext) error {
	return NewEvaluator(params).LinearTransform(ctIn, matrix, ctOut)
}
