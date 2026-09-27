package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// DropLevelSameLift drops the logical Level without changing the represented
// integer lift. It rejects transactionally when the centered source lift does
// not fit in the target Q-prefix. No coefficient division is performed.
func (eval *Evaluator) DropLevelSameLift(op0 *rlwe.Ciphertext, targetLevel int, opOut *rlwe.Ciphertext) error {
	ringQ, sourceRows, targetRows, err := eval.validateLevelDrop(op0, targetLevel, opOut)
	if err != nil {
		return err
	}
	if targetLevel < op0.Level() {
		scratch := &eval.rescaleScratch
		if scratch.coeff.N() != ringQ.N() || len(scratch.coeff.Coeffs) < sourceRows {
			eval.rescaleScratch = newFastRescaleScratch(ringQ)
			scratch = &eval.rescaleScratch
		}
		if err := scratch.validateWidth(sourceRows); err != nil {
			return err
		}
		for component := range op0.Value {
			if err := prefixToCoefficientRows(ringQ, op0.Value[component], sourceRows, op0.IsNTT, op0.IsMontgomery, scratch.coeff); err != nil {
				return fmt.Errorf("SameLift source component %d: %w", component, err)
			}
			for coefficient := 0; coefficient < ringQ.N(); coefficient++ {
				var residues [MaxQPrefixWidth]uint64
				for row := 0; row < sourceRows; row++ {
					residues[row] = scratch.coeff.Coeffs[row][coefficient]
				}
				value := reconstructQPrefix(sourceRows, residues, scratch)
				magnitude, _ := centeredQPrefix(value, scratch.modulus[sourceRows-1], scratch.half[sourceRows-1])
				if cmp192(magnitude, scratch.half[targetRows-1]) > 0 {
					return &QPrefixCapacityError{
						Level:         targetLevel,
						Component:     component,
						Bound:         uint192ToBigInt(magnitude),
						PrefixProduct: uint192ToBigInt(scratch.modulus[targetRows-1]),
					}
				}
			}
		}
	}
	commitLevelPrefix(op0, targetLevel, targetRows, opOut, eval.Parameters.N())
	return nil
}

// DropLevelCanonical explicitly forgets the source logical modulus and keeps
// the canonical target logical class. Its maintained target residues are the
// source residues modulo q0..q_min(targetLevel,3); unlike SameLift, it does not
// claim that the old centered integer lift fits in the smaller prefix.
func (eval *Evaluator) DropLevelCanonical(op0 *rlwe.Ciphertext, targetLevel int, opOut *rlwe.Ciphertext) error {
	_, _, targetRows, err := eval.validateLevelDrop(op0, targetLevel, opOut)
	if err != nil {
		return err
	}
	commitLevelPrefix(op0, targetLevel, targetRows, opOut, eval.Parameters.N())
	return nil
}

func (eval *Evaluator) validateLevelDrop(op0 *rlwe.Ciphertext, targetLevel int, opOut *rlwe.Ciphertext) (ringQ *ring.Ring, sourceRows, targetRows int, err error) {
	if eval == nil || op0 == nil || opOut == nil {
		return nil, 0, 0, errors.New("Fast Level transition evaluator and ciphertexts cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return nil, 0, 0, errors.New("Fast Level transition metadata cannot be nil")
	}
	if targetLevel < 0 || targetLevel > op0.Level() {
		return nil, 0, 0, fmt.Errorf("target Level %d must be between 0 and source Level %d", targetLevel, op0.Level())
	}
	ringQ = eval.Parameters.RingQ()
	if ringQ == nil || op0.Level() > ringQ.Level() {
		return nil, 0, 0, errors.New("Fast Level transition source exceeds configured ring Level")
	}
	if op0.N() != eval.Parameters.N() || opOut.N() != eval.Parameters.N() {
		return nil, 0, 0, errors.New("Fast Level transition dimensions do not match parameters")
	}
	if len(op0.Value) == 0 {
		return nil, 0, 0, errors.New("Fast Level transition ciphertext has no components")
	}
	sourceRows, err = QPrefixWidth(op0.Level())
	if err != nil {
		return nil, 0, 0, err
	}
	targetRows, err = QPrefixWidth(targetLevel)
	if err != nil {
		return nil, 0, 0, err
	}
	for component := range op0.Value {
		if err = validatePrefixRows(ringQ, op0.Level(), sourceRows, op0.Value[component]); err != nil {
			return nil, 0, 0, fmt.Errorf("Fast Level transition component %d: %w", component, err)
		}
	}
	return ringQ, sourceRows, targetRows, nil
}

func commitLevelPrefix(op0 *rlwe.Ciphertext, targetLevel, targetRows int, opOut *rlwe.Ciphertext, n int) {
	if opOut != op0 {
		Resize(opOut, op0.Degree(), targetLevel, n)
		for component := range op0.Value {
			copyPrefixRowsUnchecked(targetRows, op0.Value[component], opOut.Value[component])
		}
	} else {
		Resize(opOut, op0.Degree(), targetLevel, n)
	}
	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale
}
