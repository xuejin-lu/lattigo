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
// without performing Trace. The target Q-prefix is materialized from the
// canonical centered q0 representative.
func (eval *FastEvaluator) modUpBasis(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil {
		return nil, errors.New("Fast Bootstrap evaluator cannot be nil")
	}
	if eval.FastCKKS == nil {
		return nil, errors.New("Fast CKKS evaluator cannot be nil")
	}
	return eval.modUpBasisAtLevel(ct, eval.Parameters.BootstrappingParameters.MaxLevel())
}

func (eval *FastEvaluator) modUpBasisAtLevel(ct *rlwe.Ciphertext, targetLevel int) (*rlwe.Ciphertext, error) {
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
	if targetLevel < 1 || targetLevel > maxLevel {
		return nil, fmt.Errorf("Fast ModUp target level %d must be in [1,%d]", targetLevel, maxLevel)
	}
	if maxLevel < 1 {
		return nil, errors.New("Fast ModUp basis requires at least q0 and q1")
	}
	rows, err := fastckks.QPrefixWidth(targetLevel)
	if err != nil {
		return nil, err
	}
	Q := ringQ.ModuliChain()
	BRCQ := ringQ.BRedConstants()
	if rows > len(Q) || rows > len(BRCQ) || rows > len(ringQ.SubRings) {
		return nil, fmt.Errorf("Fast ModUp target prefix width %d exceeds available q rows", rows)
	}
	for component := range ct.Value {
		if len(ct.Value[component].Coeffs[0]) != params.N() {
			return nil, fmt.Errorf("Fast ModUp basis component %d has invalid q0 storage", component)
		}
	}

	// All structural and domain checks precede the first destructive transform.
	fastckks.Resize(ct, ct.Degree(), targetLevel, params.N())
	for component := range ct.Value {
		ringQ.SubRings[0].INTT(ct.Value[component].Coeffs[0], ct.Value[component].Coeffs[0])
	}

	q0 := Q[0]
	for component := range ct.Value {
		coeffs := ct.Value[component].Coeffs
		for j, residue := range coeffs[0] {
			magnitude, negative := centeredModUpQ0Residue(residue, q0)
			for row := 1; row < rows; row++ {
				q := Q[row]
				tmp := ring.BRedAdd(magnitude, q, BRCQ[row])
				if negative && tmp != 0 {
					tmp = q - tmp
				}
				coeffs[row][j] = tmp
			}
		}
	}

	for component := range ct.Value {
		for row := 0; row < rows; row++ {
			ringQ.SubRings[row].NTT(ct.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
		}
	}

	if scale := (eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()) / ct.Scale.Float64(); scale > 1 {
		scalar := uint64(math.Round(scale))
		if err := eval.FastCKKS.MulIntegerQPrefixRows(ct, new(big.Int).SetUint64(scalar), rows, ct); err != nil {
			return nil, err
		}
		ct.Scale = ct.Scale.Mul(rlwe.NewScale(scale))
	}

	return ct, nil
}

// centeredModUpQ0Residue returns the magnitude and sign of the canonical
// q0-centered representative. The midpoint q0>>1 is deliberately positive.
func centeredModUpQ0Residue(residue, q0 uint64) (magnitude uint64, negative bool) {
	if residue <= q0>>1 {
		return residue, false
	}
	return q0 - residue, true
}

// ModUp completes the Fast ModUp boundary: basis raise, explicit-width Trace,
// and consistent Montgomery conversion of every authoritative Q-prefix row.
// Dense/sparse switching and the rest of Bootstrap orchestration are outside.
func (eval *FastEvaluator) ModUp(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if eval == nil {
		return nil, errors.New("Fast Bootstrap evaluator cannot be nil")
	}
	if eval.FastCKKS == nil {
		return nil, errors.New("Fast CKKS evaluator cannot be nil")
	}
	if ct == nil || ct.MetaData == nil {
		return nil, errors.New("Fast ModUp ciphertext and metadata cannot be nil")
	}
	logSlots := eval.Parameters.CoeffsToSlotsParameters.LogSlots
	logN := eval.Parameters.BootstrappingParameters.LogN()
	if logSlots < 0 || logSlots >= logN {
		return nil, fmt.Errorf("Fast ModUp Trace logSlots must be in [0,%d), got %d", logN, logSlots)
	}
	gap := 1 << (logN - logSlots - 1)
	if logSlots == 0 {
		gap <<= 1
	}
	if gap > 1 {
		targetLevel := eval.Parameters.BootstrappingParameters.MaxLevel()
		nInv := new(big.Int).SetUint64(uint64(gap))
		if nInv.ModInverse(nInv, eval.Parameters.BootstrappingParameters.RingQ().ModulusAtLevel[targetLevel]) == nil {
			return nil, fmt.Errorf("Fast ModUp Trace cannot invert gap %d modulo Q", gap)
		}
	}
	ct, err := eval.modUpBasis(ct)
	if err != nil {
		return nil, err
	}
	rows, err := fastckks.QPrefixWidth(ct.Level())
	if err != nil {
		return nil, err
	}
	if err = eval.FastCKKS.TraceQPrefixRows(ct, logSlots, rows, ct); err != nil {
		return nil, err
	}
	ringQ := eval.Parameters.BootstrappingParameters.RingQ()
	for component := range ct.Value {
		for limb := 0; limb < rows; limb++ {
			ringQ.SubRings[limb].MForm(ct.Value[component].Coeffs[limb], ct.Value[component].Coeffs[limb])
		}
	}
	ct.IsMontgomery = true
	return ct, nil
}
