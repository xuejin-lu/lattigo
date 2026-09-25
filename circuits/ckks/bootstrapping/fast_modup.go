package bootstrapping

import (
	"errors"
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

// modUpBasis raises a Level-0 Fast ciphertext to the Bootstrap modulus basis
// through the private-F canonicalization bridge without performing Trace.
// Only maintained LogicalQ rows are exported; dormant rows remain unmaterialized.
func (eval *FastEvaluator) modUpBasis(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil || eval.FastCKKS == nil {
		return nil, errors.New("Fast Bootstrap evaluator cannot be nil")
	}
	if ct == nil || ct.MetaData == nil {
		return nil, errors.New("Fast ModUp basis ciphertext and metadata cannot be nil")
	}
	params := eval.Parameters.BootstrappingParameters
	if params.RingType() != ring.Standard {
		return nil, errors.New("Fast ModUp basis requires the Standard ring")
	}
	if ct.N() != params.N() || ct.Degree() != 1 {
		return nil, errors.New("Fast ModUp basis requires a matching degree-one ciphertext")
	}
	if ct.Level() != 0 {
		return nil, errors.New("Fast ModUp basis requires a Level-0 ciphertext")
	}
	if !ct.IsNTT {
		return nil, errors.New("Fast ModUp basis requires an NTT-domain ciphertext")
	}
	if ct.IsMontgomery {
		return nil, errors.New("Fast ModUp basis does not support Montgomery ciphertexts")
	}
	if len(ct.Value) != 2 || len(ct.Value[0].Coeffs) != 1 || len(ct.Value[1].Coeffs) != 1 {
		return nil, errors.New("Fast ModUp basis requires q0-only storage")
	}
	if ct.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return nil, errors.New("Fast ModUp basis requires a positive ciphertext scale")
	}

	maxLevel := params.MaxLevel()
	if maxLevel < 1 {
		return nil, errors.New("Fast ModUp basis requires at least q0 and q1")
	}
	for component := range ct.Value {
		if len(ct.Value[component].Coeffs[0]) != params.N() {
			return nil, fmt.Errorf("Fast ModUp basis component %d has invalid q0 storage", component)
		}
	}

	compact, err := fastckks.FusedLevel0ModUpToCompactLogical(params, ct, maxLevel, fastckks.FastCiphertextDomain{IsNTT: true})
	if err != nil {
		return nil, fmt.Errorf("Fast ModUp basis fused private-F boundary: %w", err)
	}

	if scale := (eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()) / compact.Scale.Float64(); scale > 1 {
		scalar := uint64(math.Round(scale))
		if err := eval.FastCKKS.MulIntegerMaintained(compact, new(big.Int).SetUint64(scalar), compact); err != nil {
			return nil, err
		}
		compact.Scale = compact.Scale.Mul(rlwe.NewScale(scale))
	}

	return compact, nil
}

// ModUp completes the Stage-A Fast ModUp boundary: basis raise, Fast Trace,
// and q0/q1-only Montgomery conversion. Dense/sparse switching and the rest
// of Bootstrap orchestration are intentionally outside this boundary.
func (eval *FastEvaluator) ModUp(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil {
		return nil, errors.New("Fast Bootstrap evaluator cannot be nil")
	}
	ct, err := eval.modUpBasis(ct)
	if err != nil {
		return nil, err
	}
	if err = eval.FastCKKS.Trace(ct, eval.Parameters.CoeffsToSlotsParameters.LogSlots, ct); err != nil {
		return nil, err
	}
	ringQ := eval.Parameters.BootstrappingParameters.RingQ()
	for component := range ct.Value {
		for limb := 0; limb < fastckks.MaintainedLimbCount(&eval.Parameters.BootstrappingParameters, ct.Level()); limb++ {
			ringQ.SubRings[limb].MForm(ct.Value[component].Coeffs[limb], ct.Value[component].Coeffs[limb])
		}
	}
	ct.IsMontgomery = true
	return ct, nil
}
