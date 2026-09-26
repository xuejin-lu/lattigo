package fast

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastStorageMulInteger multiplies each private-F residue by an integer
// scalar. It returns a new ciphertext and preserves the CKKS Scale metadata.
func FastStorageMulInteger(op *FastCiphertext, scalar *big.Int) (*FastCiphertext, error) {
	if err := validateFastStorageState(op); err != nil {
		return nil, err
	}
	if scalar == nil {
		return nil, fmt.Errorf("Fast storage integer scalar cannot be nil")
	}

	absScalar := new(big.Int).Abs(new(big.Int).Set(scalar))
	bounds := make([]*big.Int, len(op.componentBounds))
	for component, bound := range op.componentBounds {
		bounds[component] = new(big.Int).Mul(bound, absScalar)
	}
	if err := validateBoundsFit(op.basis, bounds, op.activeStorageWidth); err != nil {
		return nil, fmt.Errorf("Fast storage integer scalar capacity: %w", err)
	}

	output := newFastCiphertextWithBasis(op.params, op.Degree(), op.logicalLevel, op.activeStorageWidth, op.basis)
	output.metadata = cloneFastMetadata(op.metadata)
	output.componentBounds = bounds
	for row := 0; row < op.activeStorageWidth; row++ {
		subring, _ := op.basis.subring(row)
		modulus := new(big.Int).SetUint64(subring.Modulus)
		scalarMod := new(big.Int).Mod(new(big.Int).Set(scalar), modulus).Uint64()
		scalarMontgomery := ring.MForm(scalarMod, subring.Modulus, subring.BRedConstant)
		for component := range op.value {
			subring.MulScalarMontgomery(op.value[component].Coeffs[row], scalarMontgomery, output.value[component].Coeffs[row])
		}
	}
	if err := validateFastStorageState(output); err != nil {
		return nil, fmt.Errorf("invalid Fast storage integer scalar output: %w", err)
	}
	return output, nil
}

// FastStorageTraceNormalized computes the CKKS Trace in width-3 private-F
// storage. It accumulates the unnormalized automorphism sum before applying
// the exact integer normalization modulo each storage prime.
func FastStorageTraceNormalized(op *FastCiphertext, logN int) (*FastCiphertext, error) {
	if err := validateFastStorageState(op); err != nil {
		return nil, err
	}
	if op.activeStorageWidth != 3 {
		return nil, fmt.Errorf("Fast storage normalized Trace requires storage width 3, got %d", op.activeStorageWidth)
	}
	if op.Degree() != 1 {
		return nil, fmt.Errorf("Fast storage normalized Trace requires degree-one ciphertexts")
	}
	if op.logicalLevel < 1 {
		return nil, fmt.Errorf("Fast storage normalized Trace requires logical Level >= 1, got %d", op.logicalLevel)
	}
	if !op.metadata.IsNTT {
		return nil, fmt.Errorf("Fast storage normalized Trace requires NTT-domain ciphertexts")
	}
	if logN < 0 || logN >= op.LogN() {
		return nil, fmt.Errorf("Fast storage normalized Trace logN must be in [0, %d), got %d", op.LogN(), logN)
	}

	gap := 1 << uint(op.LogN()-logN-1)
	if logN == 0 {
		gap <<= 1
	}
	traceBounds := make([]*big.Int, len(op.componentBounds))
	gapBig := new(big.Int).SetUint64(uint64(gap))
	for component, bound := range op.componentBounds {
		traceBounds[component] = new(big.Int).Mul(bound, gapBig)
	}
	if err := validateBoundsFit(op.basis, traceBounds, 3); err != nil {
		return nil, fmt.Errorf("Fast storage normalized Trace capacity failure (requires strict 2*g*B < S3, g=%d): %w", gap, err)
	}

	var indexes [][]uint64
	if gap > 1 {
		logicalRing := op.params.RingQ()
		galEls := make([]uint64, 0, op.LogN()-logN)
		for i := logN; i < op.LogN()-1; i++ {
			galEls = append(galEls, op.params.GaloisElement(1<<uint(i)))
		}
		if logN == 0 {
			galEls = append(galEls, logicalRing.NthRoot()-1)
		}
		indexes = make([][]uint64, len(galEls))
		for i, galEl := range galEls {
			index, err := ring.AutomorphismNTTIndex(op.N(), logicalRing.NthRoot(), galEl)
			if err != nil {
				return nil, fmt.Errorf("Fast storage normalized Trace automorphism %d: %w", i, err)
			}
			indexes[i] = index
		}
	}

	output := op.CopyNew()
	if gap > 1 {
		permuted := make([]uint64, op.N())
		for _, index := range indexes {
			for component := range output.value {
				for row := 0; row < 3; row++ {
					subring, _ := output.basis.subring(row)
					for coefficient, source := range index {
						permuted[coefficient] = output.value[component].Coeffs[row][source]
					}
					subring.Add(output.value[component].Coeffs[row], permuted, output.value[component].Coeffs[row])
				}
			}
		}

		for row := 0; row < 3; row++ {
			subring, _ := output.basis.subring(row)
			inverse, ok := storageInverseMod(uint64(gap)%subring.Modulus, subring.Modulus)
			if !ok {
				return nil, fmt.Errorf("Fast storage normalized Trace cannot invert gap %d modulo F%d", gap, row)
			}
			scalarMontgomery := ring.MForm(inverse, subring.Modulus, subring.BRedConstant)
			for component := range output.value {
				subring.MulScalarMontgomery(output.value[component].Coeffs[row], scalarMontgomery, output.value[component].Coeffs[row])
			}
		}
	}

	if err := validateFastStorageState(output); err != nil {
		return nil, fmt.Errorf("invalid Fast storage normalized Trace output: %w", err)
	}
	return output, nil
}
