package fast

import (
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// FastStorageRescale applies one CKKS logical-Q Rescale step to a private-F
// ciphertext, returning a new value without mutating the input. The logical
// modulus q_ell is used for both coefficient division and Scale division.
func FastStorageRescale(op *FastCiphertext) (*FastCiphertext, error) {
	if err := validateFastStorageState(op); err != nil {
		return nil, err
	}
	if op.logicalLevel < 1 {
		return nil, fmt.Errorf("Fast storage Rescale requires logical Level >= 1, got %d", op.logicalLevel)
	}

	logicalQ := op.params.Q()
	q := logicalQ[op.logicalLevel]
	if q == 0 || q&1 == 0 {
		return nil, fmt.Errorf("logical Rescale divisor q_%d must be a positive odd modulus", op.logicalLevel)
	}

	output := newFastCiphertextWithBasis(op.params, op.Degree(), op.logicalLevel-1, op.activeStorageWidth, op.basis)
	output.metadata = cloneFastMetadata(op.metadata)
	output.metadata.Scale = cloneFastScale(op.metadata.Scale.Div(rlwe.NewScale(q)))
	output.metadata.IsMontgomery = false
	output.componentBounds = make([]*big.Int, len(op.componentBounds))
	for component, bound := range op.componentBounds {
		output.componentBounds[component] = rescaleComponentBound(bound, q)
	}

	coeffRows := make([][]uint64, op.activeStorageWidth)
	if op.metadata.IsNTT {
		for row := range coeffRows {
			coeffRows[row] = make([]uint64, op.N())
		}
	}
	for component := range op.value {
		for row := 0; row < op.activeStorageWidth; row++ {
			if op.metadata.IsNTT {
				subring, _ := op.basis.subring(row)
				subring.INTT(op.value[component].Coeffs[row], coeffRows[row])
			} else {
				coeffRows[row] = op.value[component].Coeffs[row]
			}
		}

		for coefficient := 0; coefficient < op.N(); coefficient++ {
			var residues [3]uint64
			for row := 0; row < op.activeStorageWidth; row++ {
				residues[row] = coeffRows[row][coefficient]
			}
			lift, err := op.basis.decodeFixed(residues, op.activeStorageWidth)
			if err != nil {
				return nil, fmt.Errorf("Fast component %d coefficient %d: %w", component, coefficient, err)
			}
			rounded, err := roundStorageIntegerByUint64(lift, q)
			if err != nil {
				return nil, fmt.Errorf("Fast component %d coefficient %d: %w", component, coefficient, err)
			}
			encoded, err := op.basis.encodeFixed(rounded, op.activeStorageWidth)
			if err != nil {
				return nil, fmt.Errorf("Fast component %d coefficient %d: %w", component, coefficient, err)
			}
			for row := 0; row < op.activeStorageWidth; row++ {
				output.value[component].Coeffs[row][coefficient] = encoded[row]
			}
		}

		if op.metadata.IsNTT {
			for row := 0; row < op.activeStorageWidth; row++ {
				subring, _ := op.basis.subring(row)
				subring.NTT(output.value[component].Coeffs[row], output.value[component].Coeffs[row])
			}
		}
	}

	if err := validateFastStorageState(output); err != nil {
		return nil, fmt.Errorf("invalid Fast storage Rescale output: %w", err)
	}
	return output, nil
}

// roundStorageIntegerByUint64 computes signed nearest-integer division for an
// odd divisor. The fixed 192-bit long division keeps per-coefficient Rescale
// free of arbitrary-precision allocations.
func roundStorageIntegerByUint64(value storageInteger, divisor uint64) (storageInteger, error) {
	if divisor == 0 {
		return storageInteger{}, fmt.Errorf("Rescale divisor cannot be zero")
	}
	if divisor&1 == 0 {
		return storageInteger{}, fmt.Errorf("nearest-integer storage division requires an odd divisor")
	}

	quotientHi, remainder := bits.Div64(0, value.magnitude.hi, divisor)
	quotientMid, remainder := bits.Div64(remainder, value.magnitude.mid, divisor)
	quotientLo, remainder := bits.Div64(remainder, value.magnitude.lo, divisor)
	quotient := storageUint192{lo: quotientLo, mid: quotientMid, hi: quotientHi}
	if remainder > divisor/2 {
		lo, carry := bits.Add64(quotient.lo, 1, 0)
		mid, carry := bits.Add64(quotient.mid, 0, carry)
		hi, overflow := bits.Add64(quotient.hi, 0, carry)
		if overflow != 0 {
			return storageInteger{}, fmt.Errorf("rounded storage quotient exceeds 192-bit representation")
		}
		quotient = storageUint192{lo: lo, mid: mid, hi: hi}
	}

	negative := value.negative && storageCmp(quotient, storageUint192{}) != 0
	return storageInteger{magnitude: quotient, negative: negative}, nil
}

func rescaleComponentBound(bound *big.Int, logicalQ uint64) *big.Int {
	numerator := new(big.Int).Set(bound)
	numerator.Add(numerator, new(big.Int).SetUint64((logicalQ-1)/2))
	return numerator.Quo(numerator, new(big.Int).SetUint64(logicalQ))
}
