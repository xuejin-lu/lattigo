package fast

import (
	"errors"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
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

// FastAddQPrefixRows adds exactly rows explicit Q-prefix limbs. Rows above the
// request are neither read nor written.
func FastAddQPrefixRows(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, rows int) error {
	return fastAddSubRows(ringQ, op0, op1, opOut, false, rows)
}

// FastSub subtracts op1 from op0 using only their authoritative q0 and q1
// limbs. Limbs q2 and above are deliberately neither read nor written.
func FastSub(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext) error {
	return fastAddSub(ringQ, op0, op1, opOut, true)
}

// FastSubQPrefixRows subtracts exactly rows explicit Q-prefix limbs. Rows
// above the request are neither read nor written.
func FastSubQPrefixRows(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, rows int) error {
	return fastAddSubRows(ringQ, op0, op1, opOut, true, rows)
}

func fastAddSub(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, sub bool) error {
	if ringQ == nil || op0 == nil || op1 == nil || opOut == nil {
		return fastAddSubRows(ringQ, op0, op1, opOut, sub, 0)
	}
	level := utils.Min(op0.Level(), op1.Level())
	level = utils.Min(level, opOut.Level())
	return fastAddSubRows(ringQ, op0, op1, opOut, sub, maintainedLimbCountForRingAtLevel(ringQ, level))
}

func fastAddSubRows(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, sub bool, rows int) error {
	return fastcore.NewAddSubWorkspace().ApplyRows(ringQ, op0, op1, opOut, rows, sub)
}

// AddQPrefixRows adds two ciphertexts using exactly rows explicit Q-prefix
// limbs. It is intended for circuits whose input producer has established a
// wider authority than the legacy Add wrapper.
func (eval *Evaluator) AddQPrefixRows(op0, op1, opOut *rlwe.Ciphertext, rows int) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if eval.addSubCore == nil {
		return errors.New("Fast Add/Sub core is not initialized")
	}
	return eval.addSubCore.ApplyRows(eval.Parameters.RingQ(), op0, op1, opOut, rows, false)
}

// SubQPrefixRows subtracts two ciphertexts using exactly rows explicit
// Q-prefix limbs. It is intended for circuits whose input producer has
// established a wider authority than the legacy Sub wrapper.
func (eval *Evaluator) SubQPrefixRows(op0, op1, opOut *rlwe.Ciphertext, rows int) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if eval.addSubCore == nil {
		return errors.New("Fast Add/Sub core is not initialized")
	}
	return eval.addSubCore.ApplyRows(eval.Parameters.RingQ(), op0, op1, opOut, rows, true)
}

func copyQ01(src, dst ring.Poly) {
	copy(dst.Coeffs[0], src.Coeffs[0])
	copy(dst.Coeffs[1], src.Coeffs[1])
}

func copyMaintained(ringQ *ring.Ring, src, dst ring.Poly) {
	level := minPolyLevel(src, dst)
	copyPrefixRowsUnchecked(maintainedLimbCountForRingAtLevel(ringQ, level), src, dst)
}
