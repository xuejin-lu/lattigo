package bootstrapping

import (
	"errors"
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// modUpBasis raises a Level-0 Fast ciphertext to the Bootstrap modulus basis
// without performing Trace. Only the maintained q0/q1 residues are computed;
// higher rows are restored structurally and remain dormant.
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

	ringQ := params.RingQ()
	maxLevel := params.MaxLevel()
	if maxLevel < 1 {
		return nil, errors.New("Fast ModUp basis requires at least q0 and q1")
	}
	ringQ0 := ringQ.AtLevel(0)
	for component := range ct.Value {
		if len(ct.Value[component].Coeffs[0]) != params.N() {
			return nil, fmt.Errorf("Fast ModUp basis component %d has invalid q0 storage", component)
		}
		ringQ0.SubRings[0].INTT(ct.Value[component].Coeffs[0], ct.Value[component].Coeffs[0])
	}

	for component := range ct.Value {
		if err := restoreFastModUpLevel(&ct.Value[component], maxLevel, params.N()); err != nil {
			return nil, fmt.Errorf("Fast ModUp basis component %d: %w", component, err)
		}
	}

	Q := ringQ.ModuliChain()
	q0 := Q[0]
	BRCQ := ringQ.BRedConstants()
	q1 := Q[1]
	for component := range ct.Value {
		coeffs := ct.Value[component].Coeffs
		for j, coeff := range coeffs[0] {
			pos, neg := uint64(1), uint64(0)
			if coeff >= q0>>1 {
				coeff = q0 - coeff
				pos, neg = 0, 1
			}
			tmp := ring.BRedAdd(coeff, q1, BRCQ[1])
			coeffs[1][j] = tmp*pos + (q1-tmp)*neg
		}
	}

	for component := range ct.Value {
		ringQ.SubRings[0].NTT(ct.Value[component].Coeffs[0], ct.Value[component].Coeffs[0])
		ringQ.SubRings[1].NTT(ct.Value[component].Coeffs[1], ct.Value[component].Coeffs[1])
	}

	if scale := (eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()) / ct.Scale.Float64(); scale > 1 {
		scalar := uint64(math.Round(scale))
		if err := eval.FastCKKS.MulIntegerMaintained(ct, new(big.Int).SetUint64(scalar), ct); err != nil {
			return nil, err
		}
		ct.Scale = ct.Scale.Mul(rlwe.NewScale(scale))
	}

	return ct, nil
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
		for limb := 0; limb < 2; limb++ {
			ringQ.SubRings[limb].MForm(ct.Value[component].Coeffs[limb], ct.Value[component].Coeffs[limb])
		}
	}
	ct.IsMontgomery = true
	return ct, nil
}

// restoreFastModUpLevel restores structural rows without allocating rows that
// were retained when ScaleDown shortened the coefficient-row slice.
func restoreFastModUpLevel(poly *ring.Poly, level, N int) error {
	if poly == nil || len(poly.Coeffs) == 0 || len(poly.Coeffs[0]) != N {
		return errors.New("invalid polynomial storage")
	}
	if cap(poly.Coeffs) < level+1 {
		coeffs := make([][]uint64, level+1)
		copy(coeffs, poly.Coeffs)
		poly.Coeffs = coeffs
	} else {
		poly.Coeffs = poly.Coeffs[:level+1]
	}
	// Fast ModUp only restores the maintained q0/q1 rows. Higher logical
	// rows remain nil so the bootstrap basis raise does not materialize the
	// dormant full-RNS storage.
	for i := 1; i <= level && i < 2; i++ {
		if len(poly.Coeffs[i]) != N {
			poly.Coeffs[i] = make([]uint64, N)
		}
	}
	return nil
}
