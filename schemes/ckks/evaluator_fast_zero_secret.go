package ckks

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
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

// mulRelinFastCKKSZeroSecret accepts only fully materialized, NTT-domain,
// non-Montgomery degree-one ciphertexts whose active c1 rows are all zero.
// Unsupported Fast ciphertext pairs fail closed rather than falling through
// to native full-Q KeySwitch/GadgetProduct. The normal CKKS product kernel
// computes the accepted c0 product into scratch without accessing a key.
func (eval Evaluator) mulRelinFastCKKSZeroSecret(op0, op1, opOut *rlwe.Ciphertext) (err error) {
	if op0 == nil || op1 == nil || opOut == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: ciphertexts cannot be nil")
	}
	if op0.Degree() != 1 || op1.Degree() != 1 {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: ciphertext inputs must both have degree 1 (got %d and %d)", op0.Degree(), op1.Degree())
	}
	if len(opOut.Value) == 0 {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: output has no polynomial storage")
	}
	if op0.MetaData == nil || op1.MetaData == nil || opOut.MetaData == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: metadata cannot be nil")
	}
	if op0.N() != eval.GetParameters().N() || op1.N() != eval.GetParameters().N() || opOut.N() != eval.GetParameters().N() {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: ciphertext ring degree does not match parameters")
	}

	ringQ := eval.GetParameters().RingQ()
	maxLevel := ringQ.Level()
	if op0.Level() < 0 || op1.Level() < 0 {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: inputs must have at least one active Q row")
	}
	if err := validateFastCKKSZeroSecretMulRows("op0", op0, op0.Level(), maxLevel, eval.GetParameters().N(), 2); err != nil {
		return err
	}
	if err := validateFastCKKSZeroSecretMulRows("op1", op1, op1.Level(), maxLevel, eval.GetParameters().N(), 2); err != nil {
		return err
	}
	if opOut.Level() < 0 {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: output has no active Q rows")
	}
	level := min(op0.Level(), op1.Level(), opOut.Level())
	outputComponents := min(len(opOut.Value), 2)
	if err := validateFastCKKSZeroSecretMulRows("output", opOut, level, maxLevel, eval.GetParameters().N(), outputComponents); err != nil {
		return err
	}

	for q := 0; q <= op0.Level(); q++ {
		for _, coefficient := range op0.Value[1].Coeffs[q] {
			if coefficient != 0 {
				return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: op0 has nonzero c1 at q%d", q)
			}
		}
	}
	for q := 0; q <= op1.Level(); q++ {
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

	// Validate every condition checked by InitOutputBinaryOp before allowing it
	// to update output metadata, so rejected inputs leave aliased outputs intact.
	if op0.IsNTT != op1.IsNTT || op0.IsNTT != eval.GetRLWEParameters().NTTFlag() {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: input NTT metadata must match the configured domain")
	}
	if op0.IsBatched != op1.IsBatched {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: input batching metadata must match")
	}
	if _, _, err := eval.InitOutputBinaryOp(op0.El(), op1.El(), 2, opOut.El()); err != nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %w", err)
	}

	// Compute into a distinct scratch ciphertext first, preserving valid in-place
	// aliasing with either input. Since both c1 components are zero, the product
	// is (a0*b0, 0, 0); relin=false never requests a key or calls GadgetProduct.
	tmpProduct := eval.pool.GetBuffCt(2, level)
	defer eval.pool.RecycleBuffCt(tmpProduct)
	tmpProduct.MetaData = &rlwe.MetaData{}
	tmpProduct.IsNTT = true
	if err := eval.mulRelin(op0, op1.El(), false, tmpProduct); err != nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %w", err)
	}

	opOut.Resize(1, level)
	opOut.Value[0].CopyLvl(level, tmpProduct.Value[0])
	opOut.Value[1].Zero()
	opOut.Scale = op0.Scale.Mul(op1.Scale)
	opOut.IsNTT = true
	opOut.IsMontgomery = false
	return nil
}

func validateFastCKKSZeroSecretMulRows(name string, ct *rlwe.Ciphertext, requiredLevel, maxLevel, n, components int) error {
	if ct == nil || len(ct.Value) < components || len(ct.Value) == 0 {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s has insufficient polynomial components", name)
	}
	level := ct.Level()
	if level < requiredLevel {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s is missing active Q rows through level %d", name, requiredLevel)
	}
	if level > maxLevel {
		return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s level %d exceeds configured Q level %d", name, level, maxLevel)
	}
	for component := 0; component < components; component++ {
		poly := ct.Value[component]
		if len(poly.Coeffs) <= requiredLevel {
			return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s component %d lacks active Q row q%d", name, component, requiredLevel)
		}
		for q := 0; q <= requiredLevel; q++ {
			if len(poly.Coeffs[q]) != n {
				return fmt.Errorf("cannot Fast CKKS zero-secret MulRelin: %s component %d q%d row is not fully materialized", name, component, q)
			}
		}
	}
	return nil
}
