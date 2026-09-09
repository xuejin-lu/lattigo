package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// FastAdd adds two ciphertexts using only their authoritative q0 and q1
// limbs. Limbs q2 and above are deliberately neither read nor written.
//
// The operands must have identical scales and q0/q1 representation domains.
// This is the Fast primitive's explicit precondition; scale alignment through
// scalar multiplication remains part of the standard evaluator and is not
// performed here.
func FastAdd(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext) error {
	return fastAddSub(ringQ, op0, op1, opOut, false)
}

// FastSub subtracts op1 from op0 using only their authoritative q0 and q1
// limbs. Limbs q2 and above are deliberately neither read nor written.
func FastSub(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext) error {
	return fastAddSub(ringQ, op0, op1, opOut, true)
}

func fastAddSub(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, sub bool) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if op0 == nil || op1 == nil || opOut == nil {
		return errors.New("op0, op1 and opOut cannot be nil")
	}
	if op0.MetaData == nil || op1.MetaData == nil || opOut.MetaData == nil {
		return errors.New("op0, op1 and opOut metadata cannot be nil")
	}
	if ringQ.Level() < 1 || op0.Level() < 1 || op1.Level() < 1 || opOut.Level() < 1 {
		return errors.New("FastAdd/FastSub requires q0 and q1")
	}
	if !op0.Scale.Equal(op1.Scale) {
		return errors.New("FastAdd/FastSub requires equal operand scales")
	}
	if op0.IsNTT != op1.IsNTT {
		return errors.New("FastAdd/FastSub requires equal IsNTT domains")
	}
	if op0.IsMontgomery != op1.IsMontgomery {
		return errors.New("FastAdd/FastSub requires equal Montgomery representations")
	}

	level := utils.Min(op0.Level(), op1.Level())
	level = utils.Min(level, opOut.Level())
	if level < 1 {
		return errors.New("FastAdd/FastSub output level must contain q0 and q1")
	}
	maxDegree := utils.Max(op0.Degree(), op1.Degree())
	minDegree := utils.Min(op0.Degree(), op1.Degree())

	validate := func(name string, ct *rlwe.Ciphertext) error {
		if ct.N() != opOut.N() || ct.N() != ringQ.N() {
			return fmt.Errorf("%s dimension does not match output", name)
		}
		for i := 0; i <= minDegree; i++ {
			if len(ct.Value[i].Coeffs) <= 1 || len(ct.Value[i].Coeffs[0]) != opOut.N() || len(ct.Value[i].Coeffs[1]) != opOut.N() {
				return fmt.Errorf("%s q0/q1 coefficient lengths are invalid", name)
			}
		}
		return nil
	}
	if err := validate("op0", op0); err != nil {
		return err
	}
	if err := validate("op1", op1); err != nil {
		return err
	}
	if opOut.N() == 0 || opOut.N() != ringQ.N() || len(opOut.Value) < maxDegree+1 {
		return errors.New("output ciphertext has invalid storage")
	}

	// Resize only changes level/degree metadata and allocation. It does not
	// cause dormant limbs to be read; the arithmetic below explicitly touches
	// q0 and q1 only.
	Resize(opOut, maxDegree, level, ringQ.N())
	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale
	opOut.IsNTT = op0.IsNTT
	opOut.IsMontgomery = op0.IsMontgomery
	opOut.IsBatched = op0.IsBatched
	opOut.LogDimensions.Rows = utils.Max(op0.LogDimensions.Rows, op1.LogDimensions.Rows)
	opOut.LogDimensions.Cols = utils.Max(op0.LogDimensions.Cols, op1.LogDimensions.Cols)

	for i := 0; i <= minDegree; i++ {
		for _, limb := range []int{0, 1} {
			dst := opOut.Value[i].Coeffs[limb]
			if sub {
				ringQ.SubRings[limb].Sub(op0.Value[i].Coeffs[limb], op1.Value[i].Coeffs[limb], dst)
			} else {
				ringQ.SubRings[limb].Add(op0.Value[i].Coeffs[limb], op1.Value[i].Coeffs[limb], dst)
			}
		}
	}

	// Preserve standard degree semantics without touching dormant limbs.
	if op0.Degree() > minDegree && opOut != op0 {
		for i := minDegree + 1; i <= op0.Degree(); i++ {
			copyQ01(op0.Value[i], opOut.Value[i])
		}
	} else if op1.Degree() > minDegree && opOut != op1 {
		for i := minDegree + 1; i <= op1.Degree(); i++ {
			if sub {
				negQ01(ringQ, op1.Value[i], opOut.Value[i])
			} else {
				copyQ01(op1.Value[i], opOut.Value[i])
			}
		}
	}

	return nil
}

func copyQ01(src, dst ring.Poly) {
	copy(dst.Coeffs[0], src.Coeffs[0])
	copy(dst.Coeffs[1], src.Coeffs[1])
}

func negQ01(ringQ *ring.Ring, src, dst ring.Poly) {
	ringQ.SubRings[0].Neg(src.Coeffs[0], dst.Coeffs[0])
	ringQ.SubRings[1].Neg(src.Coeffs[1], dst.Coeffs[1])
}
