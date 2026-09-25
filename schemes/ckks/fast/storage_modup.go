package fast

import (
	"fmt"
	"math/big"
)

// FastStorageModUpLevel0 canonicalizes a width-3 private-F ciphertext at
// logical Level 0 modulo q0, then raises only its logical Level metadata.
// The input is never modified and the output remains in the same private-F
// basis and coefficient/NTT domain.
func FastStorageModUpLevel0(op *FastCiphertext, targetLevel int) (*FastCiphertext, error) {
	if err := validateFastStorageState(op); err != nil {
		return nil, err
	}
	if op.logicalLevel != 0 {
		return nil, fmt.Errorf("Fast storage ModUp requires logical Level 0, got %d", op.logicalLevel)
	}
	if op.activeStorageWidth != 3 {
		return nil, fmt.Errorf("Fast storage ModUp requires storage width 3, got %d", op.activeStorageWidth)
	}
	if targetLevel < 1 || targetLevel > op.params.MaxLevel() {
		return nil, fmt.Errorf("target logical level %d is outside [1,%d]", targetLevel, op.params.MaxLevel())
	}

	q0 := op.params.Q()[0]
	if q0 == 0 || q0&1 == 0 {
		return nil, fmt.Errorf("logical q0 must be a positive odd modulus")
	}

	output := newFastCiphertextWithBasis(op.params, op.Degree(), targetLevel, 3, op.basis)
	output.metadata = cloneFastMetadata(op.metadata)
	output.componentBounds = make([]*big.Int, len(op.value))
	coeffRows := make([][]uint64, 3)
	if op.metadata.IsNTT {
		for row := range coeffRows {
			coeffRows[row] = make([]uint64, op.N())
		}
	}

	for component := range op.value {
		var maxMagnitude storageUint192
		for row := 0; row < 3; row++ {
			if op.metadata.IsNTT {
				subring, _ := op.basis.subring(row)
				subring.INTT(op.value[component].Coeffs[row], coeffRows[row])
			} else {
				coeffRows[row] = op.value[component].Coeffs[row]
			}
		}

		for coefficient := 0; coefficient < op.N(); coefficient++ {
			var residues [3]uint64
			for row := range residues {
				residues[row] = coeffRows[row][coefficient]
			}
			lift, err := op.basis.decodeFixed(residues, 3)
			if err != nil {
				return nil, fmt.Errorf("Fast component %d coefficient %d: %w", component, coefficient, err)
			}
			canonical := centerStorageIntegerModuloUint64(lift, q0)
			if storageCmp(canonical.magnitude, maxMagnitude) > 0 {
				maxMagnitude = canonical.magnitude
			}
			encoded, err := encodeCenteredStorageValue(op.basis, canonical, 3)
			if err != nil {
				return nil, fmt.Errorf("Fast component %d coefficient %d: %w", component, coefficient, err)
			}
			for row := 0; row < 3; row++ {
				output.value[component].Coeffs[row][coefficient] = encoded[row]
			}
		}

		output.componentBounds[component] = storageToBig(maxMagnitude)
		if op.metadata.IsNTT {
			for row := 0; row < 3; row++ {
				subring, _ := op.basis.subring(row)
				subring.NTT(output.value[component].Coeffs[row], output.value[component].Coeffs[row])
			}
		}
	}

	if err := validateFastStorageState(output); err != nil {
		return nil, fmt.Errorf("invalid Fast storage ModUp output: %w", err)
	}
	return output, nil
}

// centerStorageIntegerModuloUint64 returns Center_q(value mod q) without
// allocating arbitrary-precision integers in the coefficient loop.
func centerStorageIntegerModuloUint64(value storageInteger, modulus uint64) storageInteger {
	residue := storageMod192(value.magnitude, modulus)
	if value.negative && residue != 0 {
		residue = modulus - residue
	}
	if residue > modulus/2 {
		return storageInteger{
			magnitude: storageUint192{lo: modulus - residue},
			negative:  true,
		}
	}
	return storageInteger{magnitude: storageUint192{lo: residue}}
}
