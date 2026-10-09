package ckks

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

// addSubFastCKKSZeroSecret dispatches public ciphertext/ciphertext Add/Sub to
// the same compact Q-prefix core used by the explicit Fast evaluator.
func (eval Evaluator) addSubFastCKKSZeroSecret(op0, op1, opOut *rlwe.Ciphertext, sub bool) error {
	operation := "Add"
	if sub {
		operation = "Sub"
	}
	if op0 == nil || op1 == nil || opOut == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: ciphertexts cannot be nil", operation)
	}
	if eval.addSubCore == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: shared Fast core is unavailable", operation)
	}
	params := eval.GetParameters()
	if params.RingType() != ring.Standard {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: requires the Standard ring, got %s", operation, params.RingType())
	}
	maxLevel := params.MaxLevel()
	if op0.Level() > maxLevel || op1.Level() > maxLevel || opOut.Level() > maxLevel {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: ciphertext level exceeds configured maximum %d", operation, maxLevel)
	}
	if op0.IsNTT != op1.IsNTT || op0.IsNTT != eval.GetRLWEParameters().NTTFlag() {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: input NTT domains must match the configured domain", operation)
	}
	if op0.IsBatched != op1.IsBatched {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: input batching metadata must match", operation)
	}
	level := min(op0.Level(), op1.Level(), opOut.Level())
	if level < 1 {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: Q-prefix Add/Sub requires q0 and q1", operation)
	}
	rows, err := fastcore.QPrefixWidth(level)
	if err != nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: %w", operation, err)
	}
	if err := validateFastCKKSZeroSecretAddOperand("op0", op0, params.N(), params.RingQ(), level, rows); err != nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: %w", operation, err)
	}
	if err := validateFastCKKSZeroSecretAddOperand("op1", op1, params.N(), params.RingQ(), level, rows); err != nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: %w", operation, err)
	}
	return eval.addSubCore.ApplyRows(params.RingQ(), op0, op1, opOut, rows, sub)
}

// fastCKKSZeroSecretAddSubEligible decides whether the ordinary public
// ciphertext/ciphertext overload can preserve its existing semantics while
// using the bounded compact Fast core. Inputs that do not meet the Fast
// zero-secret contract remain on the existing full-Q CKKS implementation;
// compact inputs are separately rejected by the full-active-row preflight.
func (eval Evaluator) fastCKKSZeroSecretAddSubEligible(op0, op1 *rlwe.Ciphertext) bool {
	if op0 == nil || op1 == nil || eval.addSubCore == nil {
		return false
	}
	params := eval.GetParameters()
	if params.RingType() != ring.Standard {
		return false
	}
	maxLevel := params.MaxLevel()
	if op0.Level() < 1 || op1.Level() < 1 || op0.Level() > maxLevel || op1.Level() > maxLevel {
		return false
	}
	if !op0.Scale.Equal(op1.Scale) || op0.IsNTT != op1.IsNTT || op0.IsNTT != eval.GetRLWEParameters().NTTFlag() {
		return false
	}
	if op0.IsMontgomery != op1.IsMontgomery || op0.IsBatched != op1.IsBatched {
		return false
	}
	level := min(op0.Level(), op1.Level())
	rows, err := fastcore.QPrefixWidth(level)
	if err != nil {
		return false
	}
	if !fastCKKSZeroSecretAddBackingValid(op0, params.N(), params.RingQ(), level, rows) || !fastCKKSZeroSecretAddBackingValid(op1, params.N(), params.RingQ(), level, rows) {
		return false
	}
	if !fastCKKSZeroSecretComponentsAreZero(op0, rows) || !fastCKKSZeroSecretComponentsAreZero(op1, rows) {
		return false
	}
	return true
}

func fastCKKSZeroSecretAddBackingValid(ct *rlwe.Ciphertext, n int, ringQ *ring.Ring, level, rows int) bool {
	if ct == nil || len(ct.Value) == 0 || ct.MetaData == nil || ct.N() != n || ct.Level() < level {
		return false
	}
	for _, poly := range ct.Value {
		if err := fastcore.ValidatePrefixRows(ringQ, level, rows, poly); err != nil {
			return false
		}
	}
	return len(ct.Value) > 0
}

func fastCKKSZeroSecretComponentsAreZero(ct *rlwe.Ciphertext, rows int) bool {
	for component := 1; component < len(ct.Value); component++ {
		for row := 0; row < rows; row++ {
			for _, coefficient := range ct.Value[component].Coeffs[row] {
				if coefficient != 0 {
					return false
				}
			}
		}
	}
	return true
}

func validateFastCKKSZeroSecretAddOperand(name string, ct *rlwe.Ciphertext, n int, ringQ *ring.Ring, level, rows int) error {
	if ct == nil || len(ct.Value) == 0 {
		return fmt.Errorf("%s has no ciphertext components", name)
	}
	if ct.MetaData == nil {
		return fmt.Errorf("%s metadata cannot be nil", name)
	}
	if ct.N() != n || ct.Level() < level {
		return fmt.Errorf("%s dimensions or logical Level do not match parameters", name)
	}
	for component := range ct.Value {
		if err := fastcore.ValidatePrefixRows(ringQ, level, rows, ct.Value[component]); err != nil {
			return fmt.Errorf("%s component %d: %w", name, component, err)
		}
		for row := 0; row < rows; row++ {
			modulus := ringQ.SubRings[row].Modulus
			for coefficient, value := range ct.Value[component].Coeffs[row] {
				if value >= modulus {
					return fmt.Errorf("%s component %d q%d coefficient %d is not a canonical residue", name, component, row, coefficient)
				}
				if component > 0 && value != 0 {
					return fmt.Errorf("%s has nonzero zero-secret component %d at q%d coefficient %d", name, component, row, coefficient)
				}
			}
		}
	}
	return nil
}

// requireFastCKKSFullActiveRows prevents the still-generic scalar/plaintext
// overloads from reading dormant rows of a compact Q-prefix ciphertext.
func (eval Evaluator) requireFastCKKSFullActiveRows(operation string, ct *rlwe.Ciphertext) error {
	if ct == nil || len(ct.Value) == 0 {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: ciphertext has no polynomial storage", operation)
	}
	if ct.MetaData == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: ciphertext metadata cannot be nil", operation)
	}
	if ct.Level() < 0 || ct.Level() > eval.GetParameters().MaxLevel() {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: ciphertext Level is outside configured Q", operation)
	}
	return eval.requireFastCKKSFullActiveElementRows(operation, ct.El())
}

func (eval Evaluator) requireFastCKKSFullActiveElementRows(operation string, element *rlwe.Element[ring.Poly]) error {
	if element == nil || len(element.Value) == 0 {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: element has no polynomial storage", operation)
	}
	if element.MetaData == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: element metadata cannot be nil", operation)
	}
	if element.Level() < 0 || element.Level() > eval.GetParameters().MaxLevel() {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: element Level is outside configured Q", operation)
	}
	for component := range element.Value {
		poly := element.Value[component]
		if len(poly.Coeffs) < element.Level()+1 {
			return fmt.Errorf("cannot Fast CKKS zero-secret %s: component %d is missing active Q rows", operation, component)
		}
		for row := 0; row <= element.Level(); row++ {
			if len(poly.Coeffs[row]) != eval.GetParameters().N() {
				return fmt.Errorf("cannot Fast CKKS zero-secret %s: component %d is missing active row q%d", operation, component, row)
			}
		}
	}
	return nil
}
