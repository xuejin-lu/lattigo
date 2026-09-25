package fast

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// FusedLevel0ModUpToCompactLogical raises an ordinary Level-0 LogicalQ
// ciphertext directly to compact logical-Q rows. The private-F width-3
// boundary is semantically fused because no private-F operation observes its
// transient state here; the same canonical centered q0 lift is reduced into
// each maintained logical modulus. The width-3 capacity check uses only the
// fixed storage-prime product and does not construct private NTT subrings.
func FusedLevel0ModUpToCompactLogical(params ckks.Parameters, source *rlwe.Ciphertext, targetLevel int, target FastCiphertextDomain) (*rlwe.Ciphertext, error) {
	if source == nil {
		return nil, fmt.Errorf("source ciphertext cannot be nil")
	}
	if source.MetaData == nil {
		return nil, fmt.Errorf("source ciphertext metadata cannot be nil")
	}
	if source.IsMontgomery {
		return nil, fmt.Errorf("Montgomery LogicalQ input is not supported")
	}
	if target.IsMontgomery {
		return nil, fmt.Errorf("Montgomery LogicalQ output is not supported")
	}
	if err := validateFastCiphertextParameters(params, source.Degree(), 0, 3); err != nil {
		return nil, err
	}
	if len(source.Value) == 0 || source.N() != params.N() {
		return nil, fmt.Errorf("source ciphertext dimensions do not match CKKS parameters")
	}
	if source.Level() != 0 {
		return nil, fmt.Errorf("fused Fast ModUp requires logical Level 0, got %d", source.Level())
	}
	if targetLevel < 1 || targetLevel > params.MaxLevel() {
		return nil, fmt.Errorf("target logical level %d is outside [1,%d]", targetLevel, params.MaxLevel())
	}
	if source.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return nil, fmt.Errorf("source ciphertext scale must be positive")
	}
	for component := range source.Value {
		if len(source.Value[component].Coeffs) != 1 || len(source.Value[component].Coeffs[0]) != params.N() {
			return nil, fmt.Errorf("source component %d does not have one logical q0 row of length N", component)
		}
	}

	q0 := params.Q()[0]
	if q0 == 0 || q0&1 == 0 {
		return nil, fmt.Errorf("logical q0 must be a positive odd modulus")
	}
	if !centeredQ0FitsStorageWidth3(q0) {
		return nil, fmt.Errorf("canonical q0 lift does not fit the fixed width-3 Fast storage capacity")
	}

	maintained := maintainedLimbCount(&params, targetLevel)
	if maintained < 2 || maintained > targetLevel+1 {
		return nil, fmt.Errorf("invalid maintained logical row count %d at level %d", maintained, targetLevel)
	}

	output := NewCiphertext(&params, source.Degree(), targetLevel)
	*output.MetaData = cloneFastMetadata(*source.MetaData)
	output.IsNTT = false
	output.IsMontgomery = false

	ringQ := params.RingQ()
	half := q0 >> 1
	for component := range source.Value {
		// The compact q0 output row doubles as the only coefficient scratch.
		coefficients := output.Value[component].Coeffs[0]
		if source.IsNTT {
			ringQ.SubRings[0].INTT(source.Value[component].Coeffs[0], coefficients)
		} else {
			copy(coefficients, source.Value[component].Coeffs[0])
		}

		for coefficient, residue := range coefficients {
			if residue >= q0 {
				return nil, fmt.Errorf("source component %d coefficient %d is not canonical modulo q0", component, coefficient)
			}
			for logicalIndex := 1; logicalIndex < maintained; logicalIndex++ {
				modulus := ringQ.SubRings[logicalIndex].Modulus
				output.Value[component].Coeffs[logicalIndex][coefficient] = centeredQ0ResidueMod(residue, q0, half, modulus)
			}
		}

		if target.IsNTT {
			for logicalIndex := 0; logicalIndex < maintained; logicalIndex++ {
				ringQ.SubRings[logicalIndex].NTT(output.Value[component].Coeffs[logicalIndex], output.Value[component].Coeffs[logicalIndex])
			}
		}
	}

	output.IsNTT = target.IsNTT
	output.IsMontgomery = false
	return output, nil
}

func centeredQ0ResidueMod(residue, q0, half, modulus uint64) uint64 {
	if residue <= half {
		return residue % modulus
	}
	reduced := (q0 - residue) % modulus
	if reduced != 0 {
		return modulus - reduced
	}
	return 0
}

func centeredQ0FitsStorageWidth3(q0 uint64) bool {
	primes := fastStoragePrimes()
	product := storageMul64(storageMul64(storageUint192{lo: primes[0]}, primes[1]), primes[2])
	twiceCenteredMagnitude, overflow := storageDouble(storageUint192{lo: q0 >> 1})
	return !overflow && storageCmp(twiceCenteredMagnitude, product) < 0
}
