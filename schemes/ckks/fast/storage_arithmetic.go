package fast

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// requiredFastStorageWidth returns the smallest private-F width that proves
// centered uniqueness for every component bound.
func requiredFastStorageWidth(basis fastStorageBasis, bounds []*big.Int) (int, error) {
	if len(bounds) == 0 {
		return 0, fmt.Errorf("component bounds cannot be empty")
	}
	for component, bound := range bounds {
		if bound == nil {
			return 0, fmt.Errorf("component %d bound cannot be nil", component)
		}
		if bound.Sign() < 0 {
			return 0, fmt.Errorf("component %d bound cannot be negative", component)
		}
	}
	for width := 1; width <= 3; width++ {
		capacity := storageToBig(basis.products[width])
		fits := true
		for _, bound := range bounds {
			doubled := new(big.Int).Lsh(new(big.Int).Set(bound), 1)
			if doubled.Cmp(capacity) >= 0 {
				fits = false
				break
			}
		}
		if fits {
			return width, nil
		}
	}
	return 0, fmt.Errorf("component bounds exceed the centered capacity of all Fast storage widths")
}

// ExpandStorage returns a new ciphertext in a strictly wider private-F basis.
// It preserves the exact centered integer lift, metadata, and component bounds.
func (ct *FastCiphertext) ExpandStorage(storageWidth int) (*FastCiphertext, error) {
	if err := validateFastStorageState(ct); err != nil {
		return nil, err
	}
	if err := validateStorageWidth(storageWidth); err != nil {
		return nil, err
	}
	if storageWidth <= ct.activeStorageWidth {
		return nil, fmt.Errorf("storage expansion requires a wider basis: current=%d target=%d", ct.activeStorageWidth, storageWidth)
	}
	if err := validateBoundsFit(ct.basis, ct.componentBounds, storageWidth); err != nil {
		return nil, err
	}
	for row := 0; row < storageWidth; row++ {
		if _, err := ct.basis.subring(row); err != nil {
			return nil, err
		}
	}

	output := newFastCiphertextWithBasis(ct.params, ct.Degree(), ct.logicalLevel, storageWidth, ct.basis)
	output.metadata = cloneFastMetadata(ct.metadata)
	output.componentBounds = cloneComponentBounds(ct.componentBounds)
	coeffRows := make([][]uint64, ct.activeStorageWidth)
	for i := range coeffRows {
		coeffRows[i] = make([]uint64, ct.N())
	}

	for component := range ct.value {
		for row := 0; row < ct.activeStorageWidth; row++ {
			if ct.metadata.IsNTT {
				subring, _ := ct.basis.subring(row)
				subring.INTT(ct.value[component].Coeffs[row], coeffRows[row])
			} else {
				copy(coeffRows[row], ct.value[component].Coeffs[row])
			}
		}
		for coefficient := 0; coefficient < ct.N(); coefficient++ {
			var residues [3]uint64
			for row := 0; row < ct.activeStorageWidth; row++ {
				residues[row] = coeffRows[row][coefficient]
			}
			lift, err := ct.basis.decodeFixed(residues, ct.activeStorageWidth)
			if err != nil {
				return nil, fmt.Errorf("Fast component %d coefficient %d: %w", component, coefficient, err)
			}
			encoded, err := ct.basis.encodeFixed(lift, storageWidth)
			if err != nil {
				return nil, fmt.Errorf("Fast component %d coefficient %d: %w", component, coefficient, err)
			}
			for row := 0; row < storageWidth; row++ {
				output.value[component].Coeffs[row][coefficient] = encoded[row]
			}
		}
	}
	if ct.metadata.IsNTT {
		for component := range output.value {
			for row := 0; row < storageWidth; row++ {
				subring, _ := ct.basis.subring(row)
				subring.NTT(output.value[component].Coeffs[row], output.value[component].Coeffs[row])
			}
		}
	}
	return output, nil
}

// FastStorageAdd computes a standalone private-F ciphertext sum. It does not
// modify either input or use the logical-Q arithmetic path.
func FastStorageAdd(left, right *FastCiphertext) (*FastCiphertext, error) {
	return fastStorageAddSub(left, right, false)
}

// FastStorageSub computes a standalone private-F ciphertext difference. It
// does not modify either input or use the logical-Q arithmetic path.
func FastStorageSub(left, right *FastCiphertext) (*FastCiphertext, error) {
	return fastStorageAddSub(left, right, true)
}

func fastStorageAddSub(left, right *FastCiphertext, subtract bool) (*FastCiphertext, error) {
	if err := validateArithmeticPair(left, right); err != nil {
		return nil, err
	}
	if !left.metadata.Scale.Equal(right.metadata.Scale) {
		return nil, fmt.Errorf("Fast storage Add/Sub requires equal CKKS scales")
	}
	if !left.metadata.PlaintextMetaData.Equal(&right.metadata.PlaintextMetaData) {
		return nil, fmt.Errorf("Fast storage Add/Sub requires compatible plaintext metadata")
	}

	degree := max(left.Degree(), right.Degree())
	bounds := zeroComponentBounds(degree + 1)
	for component := range bounds {
		if component < len(left.componentBounds) {
			bounds[component].Add(bounds[component], left.componentBounds[component])
		}
		if component < len(right.componentBounds) {
			bounds[component].Add(bounds[component], right.componentBounds[component])
		}
	}
	width, err := planOperationWidth(left, right, bounds)
	if err != nil {
		return nil, err
	}
	leftRows, err := storageRowsAtWidth(left, width)
	if err != nil {
		return nil, err
	}
	rightRows, err := storageRowsAtWidth(right, width)
	if err != nil {
		return nil, err
	}

	output := newFastCiphertextWithBasis(left.params, degree, min(left.logicalLevel, right.logicalLevel), width, left.basis)
	output.metadata = cloneFastMetadata(left.metadata)
	output.componentBounds = bounds
	for component := 0; component <= degree; component++ {
		for row := 0; row < width; row++ {
			subring, _ := output.basis.subring(row)
			var leftPoly, rightPoly []uint64
			if component < len(leftRows) {
				leftPoly = leftRows[component].Coeffs[row]
			}
			if component < len(rightRows) {
				rightPoly = rightRows[component].Coeffs[row]
			}
			switch {
			case leftPoly == nil && rightPoly == nil:
				// The newly allocated output component is already zero.
			case leftPoly == nil:
				if subtract {
					subring.Sub(output.value[component].Coeffs[row], rightPoly, output.value[component].Coeffs[row])
				} else {
					copy(output.value[component].Coeffs[row], rightPoly)
				}
			case rightPoly == nil:
				copy(output.value[component].Coeffs[row], leftPoly)
			default:
				if subtract {
					subring.Sub(leftPoly, rightPoly, output.value[component].Coeffs[row])
				} else {
					subring.Add(leftPoly, rightPoly, output.value[component].Coeffs[row])
				}
			}
		}
	}
	return output, nil
}

// FastStorageMul computes raw private-F ciphertext multiplication in the NTT
// domain. It does not relinearize, rescale, or normalize the output scale.
func FastStorageMul(left, right *FastCiphertext) (*FastCiphertext, error) {
	if err := validateArithmeticPair(left, right); err != nil {
		return nil, err
	}
	if !left.metadata.IsNTT || !right.metadata.IsNTT {
		return nil, fmt.Errorf("Fast storage raw Mul requires NTT-domain operands")
	}
	if !compatiblePlaintextInterpretation(left.metadata, right.metadata) {
		return nil, fmt.Errorf("Fast storage Mul requires compatible plaintext metadata")
	}

	degree := left.Degree() + right.Degree()
	bounds := zeroComponentBounds(degree + 1)
	for outputComponent := range bounds {
		first := max(0, outputComponent-right.Degree())
		last := min(left.Degree(), outputComponent)
		for leftComponent := first; leftComponent <= last; leftComponent++ {
			rightComponent := outputComponent - leftComponent
			term := new(big.Int).Mul(left.componentBounds[leftComponent], right.componentBounds[rightComponent])
			bounds[outputComponent].Add(bounds[outputComponent], term)
		}
		bounds[outputComponent].Mul(bounds[outputComponent], new(big.Int).SetUint64(uint64(left.N())))
	}
	width, err := planOperationWidth(left, right, bounds)
	if err != nil {
		return nil, err
	}
	leftRows, err := storageRowsAtWidth(left, width)
	if err != nil {
		return nil, err
	}
	rightRows, err := storageRowsAtWidth(right, width)
	if err != nil {
		return nil, err
	}

	output := newFastCiphertextWithBasis(left.params, degree, min(left.logicalLevel, right.logicalLevel), width, left.basis)
	output.metadata = cloneFastMetadata(left.metadata)
	output.metadata.Scale = cloneFastScale(left.metadata.Scale.Mul(right.metadata.Scale))
	output.componentBounds = bounds
	for outputComponent := range output.value {
		first := max(0, outputComponent-right.Degree())
		last := min(left.Degree(), outputComponent)
		for leftComponent := first; leftComponent <= last; leftComponent++ {
			rightComponent := outputComponent - leftComponent
			for row := 0; row < width; row++ {
				subring, _ := output.basis.subring(row)
				subring.MulCoeffsBarrettThenAdd(
					leftRows[leftComponent].Coeffs[row],
					rightRows[rightComponent].Coeffs[row],
					output.value[outputComponent].Coeffs[row],
				)
			}
		}
	}
	return output, nil
}

func validateArithmeticPair(left, right *FastCiphertext) error {
	if err := validateFastStorageState(left); err != nil {
		return fmt.Errorf("left operand: %w", err)
	}
	if err := validateFastStorageState(right); err != nil {
		return fmt.Errorf("right operand: %w", err)
	}
	if left.N() != right.N() || !left.params.Equal(&right.params) {
		return fmt.Errorf("Fast storage operands require matching CKKS parameters and ring degree")
	}
	if left.metadata.IsNTT != right.metadata.IsNTT {
		return fmt.Errorf("Fast storage operands require matching coefficient/NTT domains")
	}
	return nil
}

func validateFastStorageState(ct *FastCiphertext) error {
	if ct == nil {
		return fmt.Errorf("Fast ciphertext cannot be nil")
	}
	if err := validateFastCiphertextParameters(ct.params, ct.Degree(), ct.logicalLevel, ct.activeStorageWidth); err != nil {
		return err
	}
	if ct.metadata.IsMontgomery {
		return fmt.Errorf("Montgomery Fast storage is not supported")
	}
	if ct.metadata.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return fmt.Errorf("Fast ciphertext scale must be positive")
	}
	if ct.basis.logN != ct.LogN() {
		return fmt.Errorf("Fast storage basis does not match ciphertext ring degree")
	}
	for row := 0; row < ct.activeStorageWidth; row++ {
		if _, err := ct.basis.subring(row); err != nil {
			return err
		}
	}
	if err := ct.validateStorageRows(); err != nil {
		return err
	}
	return ct.validateBoundInvariant()
}

func validateBoundsFit(basis fastStorageBasis, bounds []*big.Int, width int) error {
	if err := validateStorageWidth(width); err != nil {
		return err
	}
	capacity := storageToBig(basis.products[width])
	for component, bound := range bounds {
		if bound == nil {
			return fmt.Errorf("component %d bound cannot be nil", component)
		}
		if bound.Sign() < 0 {
			return fmt.Errorf("component %d bound cannot be negative", component)
		}
		doubled := new(big.Int).Lsh(new(big.Int).Set(bound), 1)
		if doubled.Cmp(capacity) >= 0 {
			return fmt.Errorf("component %d bound does not fit uniquely in storage width %d", component, width)
		}
	}
	return nil
}

func planOperationWidth(left, right *FastCiphertext, bounds []*big.Int) (int, error) {
	required, err := requiredFastStorageWidth(left.basis, bounds)
	if err != nil {
		return 0, err
	}
	return max(left.activeStorageWidth, right.activeStorageWidth, required), nil
}

func storageRowsAtWidth(ct *FastCiphertext, width int) ([]ring.Poly, error) {
	if ct.activeStorageWidth == width {
		return ct.value, nil
	}
	returnExpanded, err := ct.ExpandStorage(width)
	if err != nil {
		return nil, err
	}
	return returnExpanded.value, nil
}

func compatiblePlaintextInterpretation(left, right rlwe.MetaData) bool {
	return left.IsBatched == right.IsBatched &&
		left.IsBitReversed == right.IsBitReversed &&
		left.LogDimensions == right.LogDimensions
}
