package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

type fastBabyRotation struct {
	rotation int
	c0       ring.Poly
	c1       ring.Poly
}

type fastLinearTransformScratch struct {
	acc0, acc1     ring.Poly
	inner0, inner1 ring.Poly
	outer0, outer1 ring.Poly
	term0, term1   ring.Poly
	rot0, rot1     ring.Poly
	baby           []fastBabyRotation
	babyIndex      map[int]int
}

func newFastLinearTransformScratch(N int) fastLinearTransformScratch {
	newPoly := func() ring.Poly { return ring.NewPoly(N, 1) }
	return fastLinearTransformScratch{
		acc0:      newPoly(),
		acc1:      newPoly(),
		inner0:    newPoly(),
		inner1:    newPoly(),
		outer0:    newPoly(),
		outer1:    newPoly(),
		term0:     newPoly(),
		term1:     newPoly(),
		rot0:      newPoly(),
		rot1:      newPoly(),
		babyIndex: make(map[int]int),
	}
}

func (scratch *fastLinearTransformScratch) prepareBabyRotations(N int, rotations []int) {
	clear(scratch.babyIndex)
	for i, rotation := range rotations {
		if i == len(scratch.baby) {
			scratch.baby = append(scratch.baby, fastBabyRotation{
				c0: ring.NewPoly(N, 1),
				c1: ring.NewPoly(N, 1),
			})
		}
		scratch.baby[i].rotation = rotation
		scratch.babyIndex[rotation] = i
	}
}

// LinearTransform evaluates a single-level diagonal linear transformation
// using only authoritative q0 and q1 limbs. N1 == 0 uses the direct diagonal
// path; BSGS transformations use the encoded-diagonal convention from
// common/lintrans: baby rotations by i, multiplication by Vec[j+i], then a
// giant rotation by j. No evaluation key or QP basis is involved.
//
// Matrix LevelQ may be higher than the current ciphertext level. The matrix's
// q0/q1 encoding is reused at the current level and the matrix Scale is not
// modified. The output may alias the input; dormant limbs are not read or
// written.
func (eval *Evaluator) LinearTransform(ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation, ctOut *rlwe.Ciphertext) error {
	if ctIn != nil && ctOut != nil && ctIn.MetaData != nil && ctOut.MetaData != nil {
		// The receiver is an output buffer. Its domain flags may still be the
		// constructor defaults; all maintained output limbs are overwritten.
		ctOut.IsNTT = ctIn.IsNTT
		ctOut.IsMontgomery = ctIn.IsMontgomery
	}
	if err := eval.validateLinearTransform(ctIn, matrix, ctOut); err != nil {
		return err
	}

	ringQ := eval.Parameters.RingQ().AtLevel(utils.Min(ctIn.Level(), ctOut.Level()))
	var acc0, acc1 ring.Poly
	var err error
	if matrix.N1 == 0 {
		acc0, acc1, err = eval.linearTransformDirect(ringQ, ctIn, matrix)
	} else {
		acc0, acc1, err = eval.linearTransformBSGS(ringQ, ctIn, matrix)
	}
	if err != nil {
		return err
	}

	level := ringQ.Level()
	Resize(ctOut, 1, level, eval.Parameters.N())
	*ctOut.MetaData = *ctIn.MetaData
	ctOut.Scale = ctIn.Scale.Mul(matrix.Scale)
	copyQ01(acc0, ctOut.Value[0])
	copyQ01(acc1, ctOut.Value[1])
	return nil
}

func (eval *Evaluator) validateLinearTransform(ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation, ctOut *rlwe.Ciphertext) error {
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
	if ctIn.Level() < 1 || ctOut.Level() < 1 {
		return fmt.Errorf("Fast LinearTransform requires q0 and q1 (input level %d, output level %d)", ctIn.Level(), ctOut.Level())
	}
	if matrix.LevelQ < ctIn.Level() {
		return fmt.Errorf("Fast LinearTransform requires matrix.LevelQ >= ciphertext level: %d < %d", matrix.LevelQ, ctIn.Level())
	}
	if !ctIn.IsNTT || !ctOut.IsNTT || !ctIn.IsMontgomery || !ctOut.IsMontgomery || !matrix.IsNTT || !matrix.IsMontgomery {
		return errors.New("Fast LinearTransform requires NTT and Montgomery representation")
	}
	if ctIn.IsMontgomery != ctOut.IsMontgomery || ctIn.IsBatched != matrix.IsBatched {
		return errors.New("Fast LinearTransform requires matching representation and batching metadata")
	}
	for d := 0; d < 2; d++ {
		if len(ctIn.Value[d].Coeffs) < 2 || len(ctOut.Value[d].Coeffs) < 2 {
			return errors.New("Fast LinearTransform requires q0/q1 ciphertext storage")
		}
	}
	if len(matrix.Vec) == 0 {
		return errors.New("Fast LinearTransform requires at least one diagonal")
	}
	return nil
}

func (eval *Evaluator) linearTransformDirect(ringQ *ring.Ring, ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation) (acc0, acc1 ring.Poly, err error) {
	scratch := &eval.linearTransformScratch
	acc0, acc1 = scratch.acc0, scratch.acc1
	rot0, rot1 := scratch.rot0, scratch.rot1
	term0, term1 := scratch.term0, scratch.term1
	first := true
	slots := 1 << matrix.LogDimensions.Cols

	for diagonal, plaintext := range matrix.Vec {
		diagonal &= slots - 1
		if err = validateFastDiagonal(ringQ, plaintext.Q); err != nil {
			return ring.Poly{}, ring.Poly{}, fmt.Errorf("diagonal %d: %w", diagonal, err)
		}
		if err = eval.rotateComponents(ringQ, ctIn, rot0, rot1, diagonal); err != nil {
			return ring.Poly{}, ring.Poly{}, fmt.Errorf("diagonal %d automorphism: %w", diagonal, err)
		}
		fastPlaintextMul(ringQ, plaintext.Q, rot0, term0)
		fastPlaintextMul(ringQ, plaintext.Q, rot1, term1)
		if first {
			copyQ01(term0, acc0)
			copyQ01(term1, acc1)
			first = false
		} else {
			addQ01(ringQ, term0, acc0)
			addQ01(ringQ, term1, acc1)
		}
	}
	return acc0, acc1, nil
}

func (eval *Evaluator) linearTransformBSGS(ringQ *ring.Ring, ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation) (acc0, acc1 ring.Poly, err error) {
	scratch := &eval.linearTransformScratch
	acc0, acc1 = scratch.acc0, scratch.acc1
	inner0, inner1 := scratch.inner0, scratch.inner1
	outer0, outer1 := scratch.outer0, scratch.outer1
	term0, term1 := scratch.term0, scratch.term1
	index, _, rotN2 := matrix.BSGSIndex()
	scratch.prepareBabyRotations(eval.Parameters.N(), rotN2)
	eval.lastBSGSBabyRotations = 0
	for _, baby := range scratch.baby[:len(rotN2)] {
		if err = eval.rotateComponents(ringQ, ctIn, baby.c0, baby.c1, baby.rotation); err != nil {
			return ring.Poly{}, ring.Poly{}, fmt.Errorf("baby rotation %d: %w", baby.rotation, err)
		}
		if baby.rotation != 0 {
			eval.lastBSGSBabyRotations++
		}
	}
	slots := 1 << matrix.LogDimensions.Cols
	firstOuter := true

	for _, j := range utils.GetSortedKeys(index) {
		firstInner := true
		for _, i := range index[j] {
			key := j + i
			plaintext, ok := matrix.Vec[key]
			if !ok {
				plaintext, ok = matrix.Vec[key-slots]
			}
			if !ok {
				return ring.Poly{}, ring.Poly{}, fmt.Errorf("missing BSGS diagonal %d", key)
			}
			if err = validateFastDiagonal(ringQ, plaintext.Q); err != nil {
				return ring.Poly{}, ring.Poly{}, fmt.Errorf("diagonal %d: %w", key, err)
			}
			baby := scratch.baby[scratch.babyIndex[i]]
			fastPlaintextMul(ringQ, plaintext.Q, baby.c0, term0)
			fastPlaintextMul(ringQ, plaintext.Q, baby.c1, term1)
			if firstInner {
				copyQ01(term0, inner0)
				copyQ01(term1, inner1)
				firstInner = false
			} else {
				addQ01(ringQ, term0, inner0)
				addQ01(ringQ, term1, inner1)
			}
		}

		if j == 0 {
			copyQ01(inner0, outer0)
			copyQ01(inner1, outer1)
		} else {
			if err = eval.fastAutomorphism(ringQ, inner0, outer0, eval.Parameters.GaloisElement(j), true); err != nil {
				return ring.Poly{}, ring.Poly{}, fmt.Errorf("giant rotation %d: %w", j, err)
			}
			if err = eval.fastAutomorphism(ringQ, inner1, outer1, eval.Parameters.GaloisElement(j), true); err != nil {
				return ring.Poly{}, ring.Poly{}, fmt.Errorf("giant rotation %d: %w", j, err)
			}
		}
		if firstOuter {
			copyQ01(outer0, acc0)
			copyQ01(outer1, acc1)
			firstOuter = false
		} else {
			addQ01(ringQ, outer0, acc0)
			addQ01(ringQ, outer1, acc1)
		}
	}
	return acc0, acc1, nil
}

func (eval *Evaluator) rotateComponents(ringQ *ring.Ring, ctIn *rlwe.Ciphertext, out0, out1 ring.Poly, rotation int) error {
	if rotation == 0 {
		copyQ01(ctIn.Value[0], out0)
		copyQ01(ctIn.Value[1], out1)
		return nil
	}
	galEl := eval.Parameters.GaloisElement(rotation)
	if err := eval.fastAutomorphism(ringQ, ctIn.Value[0], out0, galEl, true); err != nil {
		return err
	}
	return eval.fastAutomorphism(ringQ, ctIn.Value[1], out1, galEl, true)
}

func addQ01(ringQ *ring.Ring, src, dst ring.Poly) {
	for limb := 0; limb < 2; limb++ {
		ringQ.SubRings[limb].Add(src.Coeffs[limb], dst.Coeffs[limb], dst.Coeffs[limb])
	}
}

// fastPlaintextMul multiplies an NTT/Montgomery plaintext diagonal by an
// NTT/Montgomery ciphertext polynomial using q0/q1 only.
func fastPlaintextMul(ringQ *ring.Ring, plaintext, ciphertext, output ring.Poly) {
	for limb := 0; limb < 2; limb++ {
		ringQ.SubRings[limb].MulCoeffsMontgomery(plaintext.Coeffs[limb], ciphertext.Coeffs[limb], output.Coeffs[limb])
	}
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
