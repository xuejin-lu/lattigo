package fast

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// FastCiphertextDomain explicitly selects the coefficient/NTT and
// ordinary/Montgomery representation at a basis-conversion boundary.
type FastCiphertextDomain struct {
	IsNTT        bool
	IsMontgomery bool
}

// ImportLevel0 converts an ordinary Standard-ring LogicalQ ciphertext at
// Level 0 into physically separate private-F rows. It uses the canonical
// centered q0 representative and rejects Montgomery input and insufficient
// active storage capacity rather than reinterpreting or wrapping rows.
func ImportLevel0(params ckks.Parameters, source *rlwe.Ciphertext, storageWidth int, target FastCiphertextDomain) (*FastCiphertext, error) {
	if source == nil {
		return nil, fmt.Errorf("source ciphertext cannot be nil")
	}
	if source.MetaData == nil {
		return nil, fmt.Errorf("source ciphertext metadata cannot be nil")
	}
	if source.IsMontgomery {
		return nil, fmt.Errorf("Montgomery LogicalQ input is not supported by Fast storage import")
	}
	if target.IsMontgomery {
		return nil, fmt.Errorf("Montgomery Fast storage output is not supported")
	}
	if err := validateFastCiphertextParameters(params, source.Degree(), 0, storageWidth); err != nil {
		return nil, err
	}
	if len(source.Value) == 0 || source.N() != params.N() {
		return nil, fmt.Errorf("source ciphertext dimensions do not match CKKS parameters")
	}
	if source.Level() != 0 {
		return nil, fmt.Errorf("Fast storage import requires logical Level 0, got %d", source.Level())
	}
	if source.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return nil, fmt.Errorf("source ciphertext scale must be positive")
	}
	for component := range source.Value {
		if len(source.Value[component].Coeffs) != 1 || len(source.Value[component].Coeffs[0]) != params.N() {
			return nil, fmt.Errorf("source component %d does not have one logical q0 row of length N", component)
		}
	}

	basis, err := newFastStorageBasis(params.LogN())
	if err != nil {
		return nil, err
	}
	fast := newFastCiphertextWithBasis(params, source.Degree(), 0, storageWidth, basis)
	fast.metadata = cloneFastMetadata(*source.MetaData)
	fast.metadata.IsNTT = target.IsNTT
	fast.metadata.IsMontgomery = false

	logicalRing := params.RingQ().AtLevel(0)
	q0 := logicalRing.SubRings[0].Modulus
	coefficients := make([]uint64, params.N())
	for component := range source.Value {
		sourceRow := source.Value[component].Coeffs[0]
		if source.IsNTT {
			logicalRing.SubRings[0].INTT(sourceRow, coefficients)
		} else {
			copy(coefficients, sourceRow)
		}
		for k, residue := range coefficients {
			if residue >= q0 {
				return nil, fmt.Errorf("source component %d coefficient %d is not canonical modulo q0", component, k)
			}
			lift := centeredStorageResidue(residue, q0)
			encoded, err := encodeCenteredStorageValue(basis, lift, storageWidth)
			if err != nil {
				return nil, fmt.Errorf("source component %d coefficient %d: %w", component, k, err)
			}
			for i := 0; i < storageWidth; i++ {
				fast.value[component].Coeffs[i][k] = encoded[i]
			}
		}
	}
	if target.IsNTT {
		for component := range fast.value {
			for i := 0; i < storageWidth; i++ {
				subring, _ := basis.subring(i)
				subring.NTT(fast.value[component].Coeffs[i], fast.value[component].Coeffs[i])
			}
		}
	}
	return fast, nil
}

func centeredStorageResidue(residue, modulus uint64) storageInteger {
	if residue > modulus/2 {
		return storageInteger{magnitude: storageUint192{lo: modulus - residue}, negative: true}
	}
	return storageInteger{magnitude: storageUint192{lo: residue}}
}

func encodeCenteredStorageValue(basis fastStorageBasis, value storageInteger, width int) ([3]uint64, error) {
	var residues [3]uint64
	if err := validateStorageWidth(width); err != nil {
		return residues, err
	}
	doubled, overflow := storageDouble(value.magnitude)
	if overflow || storageCmp(doubled, basis.products[width]) >= 0 {
		return residues, fmt.Errorf("centered lift does not fit uniquely in storage width %d", width)
	}
	return basis.encodeFixed(value, width)
}

func storageDouble(value storageUint192) (storageUint192, bool) {
	lo, carry := bits.Add64(value.lo, value.lo, 0)
	mid, carry := bits.Add64(value.mid, value.mid, carry)
	hi, overflow := bits.Add64(value.hi, value.hi, carry)
	return storageUint192{lo: lo, mid: mid, hi: hi}, overflow != 0
}

// ExportToLogical reconstructs each authoritative centered Fast lift in the
// coefficient domain, reduces it into each logical q_i row through the
// ciphertext's explicit LogicalLevel, and optionally applies logical NTT.
// Montgomery export is explicitly unsupported.
func (ct *FastCiphertext) ExportToLogical(target FastCiphertextDomain) (*rlwe.Ciphertext, error) {
	if ct == nil {
		return nil, fmt.Errorf("Fast ciphertext cannot be nil")
	}
	if target.IsMontgomery {
		return nil, fmt.Errorf("Montgomery LogicalQ output is not supported")
	}
	if ct.metadata.IsMontgomery {
		return nil, fmt.Errorf("Montgomery Fast storage input is not supported")
	}
	if err := validateFastCiphertextParameters(ct.params, ct.Degree(), ct.logicalLevel, ct.activeStorageWidth); err != nil {
		return nil, err
	}
	if err := ct.validateStorageRows(); err != nil {
		return nil, err
	}

	output := rlwe.NewCiphertext(ct.params, ct.Degree(), ct.logicalLevel)
	*output.MetaData = cloneFastMetadata(ct.metadata)
	output.IsNTT = false
	output.IsMontgomery = false

	coeffRows := make([][]uint64, ct.activeStorageWidth)
	for i := range coeffRows {
		coeffRows[i] = make([]uint64, ct.N())
	}
	for component := range ct.value {
		for i := 0; i < ct.activeStorageWidth; i++ {
			if ct.metadata.IsNTT {
				subring, _ := ct.basis.subring(i)
				subring.INTT(ct.value[component].Coeffs[i], coeffRows[i])
			} else {
				copy(coeffRows[i], ct.value[component].Coeffs[i])
			}
		}
		for k := 0; k < ct.N(); k++ {
			var storageResidues [3]uint64
			for i := 0; i < ct.activeStorageWidth; i++ {
				storageResidues[i] = coeffRows[i][k]
			}
			lift, err := ct.basis.decodeFixed(storageResidues, ct.activeStorageWidth)
			if err != nil {
				return nil, fmt.Errorf("Fast component %d coefficient %d: %w", component, k, err)
			}
			for logicalIndex := 0; logicalIndex <= ct.logicalLevel; logicalIndex++ {
				q := ct.params.RingQ().SubRings[logicalIndex].Modulus
				residue := storageMod192(lift.magnitude, q)
				if lift.negative && residue != 0 {
					residue = q - residue
				}
				output.Value[component].Coeffs[logicalIndex][k] = residue
			}
		}
	}
	if target.IsNTT {
		logicalRing := ct.params.RingQ().AtLevel(ct.logicalLevel)
		for component := range output.Value {
			for logicalIndex := 0; logicalIndex <= ct.logicalLevel; logicalIndex++ {
				logicalRing.SubRings[logicalIndex].NTT(output.Value[component].Coeffs[logicalIndex], output.Value[component].Coeffs[logicalIndex])
			}
		}
	}
	output.IsNTT = target.IsNTT
	output.IsMontgomery = false
	return output, nil
}

func (ct *FastCiphertext) validateStorageRows() error {
	if ct.activeStorageWidth < 1 || ct.activeStorageWidth > 3 {
		return fmt.Errorf("invalid active storage width %d", ct.activeStorageWidth)
	}
	if len(ct.value) == 0 {
		return fmt.Errorf("Fast ciphertext must contain at least one component")
	}
	primes := fastStoragePrimes()
	for component := range ct.value {
		if len(ct.value[component].Coeffs) != ct.activeStorageWidth {
			return fmt.Errorf("Fast component %d has %d rows, expected storage width %d", component, len(ct.value[component].Coeffs), ct.activeStorageWidth)
		}
		for i := 0; i < ct.activeStorageWidth; i++ {
			row := ct.value[component].Coeffs[i]
			if len(row) != ct.N() {
				return fmt.Errorf("Fast component %d storage row %d has length %d, expected N=%d", component, i, len(row), ct.N())
			}
			for k, residue := range row {
				if residue >= primes[i] {
					return fmt.Errorf("Fast component %d storage row %d coefficient %d is not canonical", component, i, k)
				}
			}
		}
	}
	return nil
}
