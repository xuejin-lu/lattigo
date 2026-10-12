package ckks

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func isFastCKKSScalarOperand(value rlwe.Operand) bool {
	switch value.(type) {
	case complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex:
		return true
	default:
		return false
	}
}

func (eval Evaluator) addSubFastCKKSZeroSecretScalar(operation string, op0 *rlwe.Ciphertext, scalar rlwe.Operand, opOut *rlwe.Ciphertext, sub bool) error {
	_, rows, err := eval.validateFastCKKSScalarPair(operation, op0, opOut, false)
	if err != nil {
		return err
	}
	value := scalar
	if uintValue, ok := value.(uint); ok {
		value = uint64(uintValue)
	}
	encoded := bignum.ToComplex(value, eval.GetParameters().EncodingPrecision())
	if eval.scalarCore == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret %s: shared scalar core is unavailable", operation)
	}
	return eval.scalarCore.AddSubRows(eval.GetParameters().RingQ(), op0, opOut, rows, encoded, sub)
}

func (eval Evaluator) mulScalarFastCKKSZeroSecret(op0 *rlwe.Ciphertext, scalar rlwe.Operand, opOut *rlwe.Ciphertext) error {
	_, rows, err := eval.validateFastCKKSScalarPair("Mul", op0, opOut, false)
	if err != nil {
		return err
	}
	value := scalar
	if uintValue, ok := value.(uint); ok {
		value = uint64(uintValue)
	}
	encoded := bignum.ToComplex(value, eval.GetParameters().EncodingPrecision())
	if eval.scalarCore == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret Mul: shared scalar core is unavailable")
	}
	return eval.scalarCore.MulRows(eval.GetParameters().RingQ(), op0, opOut, rows, eval.GetParameters().LevelsConsumedPerRescaling(), encoded)
}

func (eval Evaluator) mulThenAddScalarFastCKKSZeroSecret(op0 *rlwe.Ciphertext, scalar rlwe.Operand, opOut *rlwe.Ciphertext) error {
	_, rows, err := eval.validateFastCKKSScalarPair("MulThenAdd", op0, opOut, true)
	if err != nil {
		return err
	}
	value := scalar
	if uintValue, ok := value.(uint); ok {
		value = uint64(uintValue)
	}
	encoded := bignum.ToComplex(value, eval.GetParameters().EncodingPrecision())
	if eval.scalarCore == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulThenAdd: shared scalar core is unavailable")
	}
	if err := eval.scalarCore.MulThenAddRows(eval.GetParameters().RingQ(), op0, opOut, rows, eval.GetParameters().LevelsConsumedPerRescaling(), encoded); err != nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulThenAdd: %w", err)
	}
	return nil
}

func (eval Evaluator) validateFastCKKSScalarPair(operation string, op0, opOut *rlwe.Ciphertext, outputIsAccumulator bool) (level, rows int, err error) {
	fail := func(format string, args ...any) (int, int, error) {
		return 0, 0, fmt.Errorf("cannot Fast CKKS zero-secret %s: %s", operation, fmt.Sprintf(format, args...))
	}
	if eval.GetParameters().RingType() != ring.Standard {
		return fail("requires the Standard ring, got %s", eval.GetParameters().RingType())
	}
	level0, err := validateFastCKKSScalarCiphertext("input", op0, eval.GetParameters(), true)
	if err != nil {
		return fail("%v", err)
	}
	levelOut, err := validateFastCKKSScalarCiphertext("output", opOut, eval.GetParameters(), true)
	if err != nil {
		return fail("%v", err)
	}
	if outputIsAccumulator && op0.IsBatched != opOut.IsBatched {
		return fail("input/output batching metadata must match")
	}
	if op0.IsNTT != eval.GetRLWEParameters().NTTFlag() || !op0.IsNTT || opOut.IsNTT != op0.IsNTT {
		return fail("operands must use the configured NTT domain")
	}
	if op0.IsMontgomery != opOut.IsMontgomery {
		return fail("input/output Montgomery representations must match")
	}
	level = min(level0, levelOut)
	rows, err = fastcore.QPrefixWidth(level)
	if err != nil {
		return fail("%v", err)
	}
	return level, rows, nil
}

func validateFastCKKSScalarCiphertext(name string, ct *rlwe.Ciphertext, params *Parameters, allowLevelZero bool) (int, error) {
	if ct == nil {
		return 0, fmt.Errorf("%s ciphertext cannot be nil", name)
	}
	if len(ct.Value) == 0 {
		return 0, fmt.Errorf("%s ciphertext has no polynomial components", name)
	}
	if ct.MetaData == nil {
		return 0, fmt.Errorf("%s ciphertext metadata cannot be nil", name)
	}
	if len(ct.Value[0].Coeffs) == 0 || len(ct.Value[0].Coeffs[0]) == 0 {
		return 0, fmt.Errorf("%s ciphertext has empty c0 backing", name)
	}
	level := len(ct.Value[0].Coeffs) - 1
	minLevel := 1
	if allowLevelZero {
		minLevel = 0
	}
	if params == nil || params.RingQ() == nil || ct.N() != params.N() || level < minLevel || level > params.MaxLevel() {
		return 0, fmt.Errorf("%s ciphertext has invalid dimensions or logical Level", name)
	}
	rows, err := fastcore.QPrefixWidth(level)
	if err != nil {
		return 0, err
	}
	if len(ct.Value) < 1 {
		return 0, fmt.Errorf("%s ciphertext has no components", name)
	}
	for component := range ct.Value {
		poly := ct.Value[component]
		if len(poly.Coeffs) != level+1 {
			return 0, fmt.Errorf("%s component %d has inconsistent logical Q rows", name, component)
		}
		if err := fastcore.ValidatePrefixRows(params.RingQ(), level, rows, poly); err != nil {
			return 0, fmt.Errorf("%s component %d: %w", name, component, err)
		}
	}
	return level, nil
}

func (eval Evaluator) mulFastCKKSPrefix(op0, op1, opOut *rlwe.Ciphertext) error {
	if op0 == nil || op1 == nil || opOut == nil {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: ciphertexts cannot be nil")
	}
	params := eval.GetParameters()
	level0, err := validateFastCKKSScalarCiphertext("op0", op0, params, false)
	if err != nil {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: %w", err)
	}
	level1, err := validateFastCKKSScalarCiphertext("op1", op1, params, false)
	if err != nil {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: %w", err)
	}
	if op0.Degree() != 1 || op1.Degree() != 1 {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: ciphertext inputs must both have degree 1")
	}
	levelOut, err := validateFastCKKSScalarCiphertext("output", opOut, params, false)
	if err != nil {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: %w", err)
	}
	if op0.IsNTT != eval.GetRLWEParameters().NTTFlag() || !op0.IsNTT || op1.IsNTT != op0.IsNTT || opOut.IsNTT != op0.IsNTT {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: operands must use the configured NTT domain")
	}
	if op0.IsMontgomery != op1.IsMontgomery || opOut.IsMontgomery != op0.IsMontgomery || op0.IsBatched != op1.IsBatched {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: operand representations and batching metadata must match")
	}
	level := min(level0, level1, levelOut)
	rows, err := fastcore.QPrefixWidth(level)
	if err != nil {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: %w", err)
	}
	if eval.mulCore == nil {
		return fmt.Errorf("cannot Fast CKKS Q-prefix Mul: shared Fast Mul core is unavailable")
	}
	return eval.mulCore.ApplyRows(params.RingQ(), op0, op1.El(), opOut, rows, false)
}
