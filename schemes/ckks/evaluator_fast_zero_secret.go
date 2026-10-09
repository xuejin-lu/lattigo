package ckks

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

// fastCKKSZeroSecretSimulationProvider is implemented only by the Fast CKKS
// parameter type. It mirrors the capability used by core/rlwe's Encryptor
// without making the public CKKS evaluator depend on the sibling fast package.
type fastCKKSZeroSecretSimulationProvider interface {
	FastCKKSZeroSecretSimulation()
}

func (eval Evaluator) usesFastCKKSZeroSecretSimulation() bool {
	_, ok := any(eval.GetParameters()).(fastCKKSZeroSecretSimulationProvider)
	return ok
}

// mulRelinFastCKKSZeroSecret accepts only NTT-domain, non-Montgomery
// degree-one ciphertexts whose authoritative c1 rows are all zero. The shared
// Fast core consumes only the maintained Q-prefix; dormant logical rows are
// neither required nor read.
func (eval Evaluator) mulRelinFastCKKSZeroSecret(op0, op1, opOut *rlwe.Ciphertext) (err error) {
	if op0 == nil || op1 == nil || opOut == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: ciphertexts cannot be nil")
	}
	params := eval.GetParameters()
	if _, err := validateFastCKKSZeroSecretMulShape("op0", op0, params); err != nil {
		return err
	}
	if _, err := validateFastCKKSZeroSecretMulShape("op1", op1, params); err != nil {
		return err
	}
	if _, err := validateFastCKKSZeroSecretMulBacking("output", opOut, params); err != nil {
		return err
	}
	if len(opOut.Value) < 2 {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: output has insufficient polynomial components")
	}
	if op0.IsNTT != op1.IsNTT {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: input NTT metadata must match")
	}
	if op0.IsBatched != op1.IsBatched {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: input batching metadata must match")
	}
	ringQ := params.RingQ()
	maxLevel := ringQ.Level()
	level := min(op0.Level(), op1.Level(), opOut.Level())
	rows, err := fastcore.QPrefixWidth(level)
	if err != nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %w", err)
	}
	if err := validateFastCKKSZeroSecretMulRows("op0", op0, level, maxLevel, params.N(), 2, rows, ringQ); err != nil {
		return err
	}
	if err := validateFastCKKSZeroSecretMulRows("op1", op1, level, maxLevel, params.N(), 2, rows, ringQ); err != nil {
		return err
	}
	outputComponents := min(len(opOut.Value), 2)
	if err := validateFastCKKSZeroSecretMulRows("output", opOut, level, maxLevel, params.N(), outputComponents, rows, ringQ); err != nil {
		return err
	}

	for q := 0; q < rows; q++ {
		modulus := ringQ.SubRings[q].Modulus
		for coefficient, value := range op0.Value[0].Coeffs[q] {
			if value >= modulus {
				return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: op0 c0 q%d coefficient %d is not a canonical residue", q, coefficient)
			}
		}
		for _, coefficient := range op0.Value[1].Coeffs[q] {
			if coefficient != 0 {
				return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: op0 has nonzero c1 at q%d", q)
			}
		}
	}
	for q := 0; q < rows; q++ {
		modulus := ringQ.SubRings[q].Modulus
		for coefficient, value := range op1.Value[0].Coeffs[q] {
			if value >= modulus {
				return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: op1 c0 q%d coefficient %d is not a canonical residue", q, coefficient)
			}
		}
		for _, coefficient := range op1.Value[1].Coeffs[q] {
			if coefficient != 0 {
				return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: op1 has nonzero c1 at q%d", q)
			}
		}
	}

	if !op0.IsNTT || !op1.IsNTT || !eval.GetRLWEParameters().NTTFlag() {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: inputs must use the configured NTT domain")
	}
	if op0.IsMontgomery || op1.IsMontgomery {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: Montgomery-form inputs are not supported")
	}

	if op0.IsNTT != op1.IsNTT || op0.IsNTT != eval.GetRLWEParameters().NTTFlag() {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: input NTT metadata must match the configured domain")
	}
	if op0.IsBatched != op1.IsBatched {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: input batching metadata must match")
	}
	if eval.mulCore == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: shared Fast Mul core is unavailable")
	}

	// The core computes into evaluator-owned scratch before committing output,
	// so valid in-place aliases remain safe and errors leave caller storage intact.
	return eval.mulCore.ApplyRows(ringQ, op0, op1.El(), opOut, rows, true)
}

// validateFastCKKSZeroSecretMulShape checks the metadata and backing needed to
// safely inspect a degree-one ciphertext's Level and dimensions. It examines
// only rows in the fixed Q-prefix policy; dormant rows remain untouched.
func validateFastCKKSZeroSecretMulShape(name string, ct *rlwe.Ciphertext, params *Parameters) (int, error) {
	level, err := validateFastCKKSZeroSecretMulBacking(name, ct, params)
	if err != nil {
		return 0, err
	}
	if len(ct.Value) != 2 {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s must have exactly two components", name)
	}
	return level, nil
}

func validateFastCKKSZeroSecretMulBacking(name string, ct *rlwe.Ciphertext, params *Parameters) (int, error) {
	if ct == nil {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s cannot be nil", name)
	}
	if len(ct.Value) == 0 {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s has no polynomial components", name)
	}
	if ct.MetaData == nil {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s metadata cannot be nil", name)
	}
	if len(ct.Value[0].Coeffs) == 0 || len(ct.Value[0].Coeffs[0]) == 0 {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s has empty polynomial backing", name)
	}
	level := len(ct.Value[0].Coeffs) - 1
	if params == nil || params.RingQ() == nil || level < 0 || level > params.MaxLevel() {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s has an invalid logical Level", name)
	}
	rows, err := fastcore.QPrefixWidth(level)
	if err != nil {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s: %w", name, err)
	}
	for component := range ct.Value {
		poly := ct.Value[component]
		if len(poly.Coeffs) != level+1 {
			return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s component %d has inconsistent logical rows", name, component)
		}
		for row := 0; row < rows; row++ {
			if len(poly.Coeffs[row]) != params.N() {
				return 0, fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s component %d q%d backing length does not match N=%d", name, component, row, params.N())
			}
		}
	}
	return level, nil
}

func validateFastCKKSZeroSecretMulRows(name string, ct *rlwe.Ciphertext, requiredLevel, maxLevel, n, components, rows int, ringQ *ring.Ring) error {
	if ct == nil || len(ct.Value) < components || len(ct.Value) == 0 {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s has insufficient polynomial components", name)
	}
	if ct.MetaData == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s metadata cannot be nil", name)
	}
	level := ct.Level()
	if level < requiredLevel {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s is missing active Q rows through level %d", name, requiredLevel)
	}
	if level > maxLevel {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s level %d exceeds configured Q level %d", name, level, maxLevel)
	}
	if ct.N() != n {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s ring degree does not match parameters", name)
	}
	if components == 0 {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s has no output components", name)
	}
	polys := ct.Value[:components]
	if err := fastcore.ValidatePrefixRows(ringQ, requiredLevel, rows, polys...); err != nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s: %w", name, err)
	}
	return nil
}
