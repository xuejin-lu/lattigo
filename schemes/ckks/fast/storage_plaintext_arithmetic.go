package fast

import (
	"errors"
	"fmt"
	"math/big"
	"sort"

	"github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// FastStorageMulPlaintext multiplies each private-F ciphertext component by
// one private-F plaintext polynomial. It performs only active-width NTT-domain
// arithmetic and does not reconstruct coefficients in the execution path.
func FastStorageMulPlaintext(ct *FastCiphertext, pt *FastStoragePlaintext) (*FastCiphertext, error) {
	if err := validateFastStoragePlaintextOperand(ct, pt); err != nil {
		return nil, err
	}
	bounds := make([]*big.Int, len(ct.componentBounds))
	for component, bound := range ct.componentBounds {
		bounds[component] = new(big.Int).Mul(bound, pt.l1Norm)
	}
	if err := validateBoundsFit(ct.basis, bounds, ct.activeStorageWidth); err != nil {
		return nil, fmt.Errorf("private-F plaintext multiplication capacity: %w", err)
	}

	output := newFastCiphertextWithBasis(ct.params, ct.Degree(), ct.logicalLevel, ct.activeStorageWidth, ct.basis)
	output.metadata = cloneFastMetadata(ct.metadata)
	output.metadata.Scale = ct.metadata.Scale.Mul(pt.scale)
	output.componentBounds = bounds
	for row := 0; row < ct.activeStorageWidth; row++ {
		subring, _ := ct.basis.subring(row)
		for component := range ct.value {
			subring.MulCoeffsBarrett(ct.value[component].Coeffs[row], pt.value.Coeffs[row], output.value[component].Coeffs[row])
		}
	}
	if err := validateFastStorageState(output); err != nil {
		return nil, fmt.Errorf("invalid private-F plaintext multiplication output: %w", err)
	}
	return output, nil
}

func validateFastStoragePlaintextOperand(ct *FastCiphertext, pt *FastStoragePlaintext) error {
	if err := validateFastStorageState(ct); err != nil {
		return fmt.Errorf("ciphertext: %w", err)
	}
	if err := validateFastStoragePlaintext(pt); err != nil {
		return fmt.Errorf("plaintext: %w", err)
	}
	if ct.activeStorageWidth != pt.width || !ct.metadata.IsNTT || ct.metadata.IsMontgomery {
		return fmt.Errorf("private-F plaintext arithmetic requires matching-width NTT ordinary-residue storage (ciphertext=%d plaintext=%d)", ct.activeStorageWidth, pt.width)
	}
	if ct.N() != pt.params.N() || !ct.params.Equal(&pt.params) || ct.basis.logN != pt.basis.logN {
		return errors.New("private-F plaintext and ciphertext require matching CKKS parameters")
	}
	if pt.levelQ < ct.logicalLevel {
		return fmt.Errorf("plaintext LevelQ %d is below ciphertext logical Level %d", pt.levelQ, ct.logicalLevel)
	}
	if ct.metadata.IsBatched != pt.metadata.IsBatched || ct.metadata.IsBitReversed != pt.metadata.IsBitReversed || ct.metadata.LogDimensions != pt.metadata.LogDimensions {
		return errors.New("private-F plaintext and ciphertext batching metadata are incompatible")
	}
	return nil
}

// FastStorageAutomorphism applies a CKKS Galois automorphism to every NTT
// private-F row and ciphertext component without key switching. As in the
// current Fast path, this is a representation/arithmetic primitive, not a
// public-key rotation operation.
func FastStorageAutomorphism(ct *FastCiphertext, galEl uint64) (*FastCiphertext, error) {
	if err := validateFastStorageState(ct); err != nil {
		return nil, err
	}
	if !ct.metadata.IsNTT || ct.metadata.IsMontgomery {
		return nil, errors.New("private-F automorphism requires NTT ordinary-residue storage")
	}
	index, err := ring.AutomorphismNTTIndex(ct.N(), ct.params.RingQ().NthRoot(), galEl)
	if err != nil {
		return nil, fmt.Errorf("private-F automorphism index: %w", err)
	}
	output := ct.CopyNew()
	for component := range ct.value {
		for row := 0; row < ct.activeStorageWidth; row++ {
			permuteStorageNTT(ct.value[component].Coeffs[row], output.value[component].Coeffs[row], index)
		}
	}
	if err := validateFastStorageState(output); err != nil {
		return nil, fmt.Errorf("invalid private-F automorphism output: %w", err)
	}
	return output, nil
}

func permuteStorageNTT(input, output []uint64, index []uint64) {
	for coefficient, source := range index {
		output[coefficient] = input[source]
	}
}

// FastStorageLinearTransformation contains immutable private-F mirrors and
// the exact direct/BSGS schedule of one encoded LogicalQ transformation.
type FastStorageLinearTransformation struct {
	params    ckks.Parameters
	width     int
	levelQ    int
	scale     rlwe.Scale
	metadata  rlwe.PlaintextMetaData
	n1        int
	cols      int
	diagonals map[int]*FastStoragePlaintext
	index     map[int][]int
	sumL1Norm *big.Int
	maxAbs    *big.Int
	maxL1Norm *big.Int
}

// NewFastStorageLinearTransformation mirrors every encoded Q diagonal and
// retains the corresponding direct/BSGS execution schedule.
func NewFastStorageLinearTransformation(params ckks.Parameters, matrix lintrans.LinearTransformation) (*FastStorageLinearTransformation, error) {
	return NewFastStorageLinearTransformationWithWidth(params, matrix, 3)
}

// NewFastStorageLinearTransformationWithWidth mirrors every encoded Q
// diagonal into the requested width-2 or width-3 private-F basis, preserving
// the exact direct/BSGS schedule.
func NewFastStorageLinearTransformationWithWidth(params ckks.Parameters, matrix lintrans.LinearTransformation, storageWidth int) (*FastStorageLinearTransformation, error) {
	if storageWidth != 2 && storageWidth != 3 {
		return nil, fmt.Errorf("private-F LinearTransformation width must be 2 or 3, got %d", storageWidth)
	}
	if matrix.MetaData == nil {
		return nil, errors.New("LogicalQ linear transformation metadata cannot be nil")
	}
	if params.N() == 0 || params.RingType() != ring.Standard || matrix.LevelQ < 0 || matrix.LevelQ > params.MaxLevel() {
		return nil, errors.New("linear transformation has invalid Standard CKKS parameters or LevelQ")
	}
	if !matrix.IsNTT || !matrix.IsMontgomery {
		return nil, errors.New("LogicalQ linear transformation diagonals must be NTT and Montgomery encoded")
	}
	if matrix.Scale.Cmp(rlwe.NewScale(0)) != 1 || len(matrix.Vec) == 0 {
		return nil, errors.New("linear transformation requires a positive Scale and at least one diagonal")
	}
	if matrix.N1 < 0 {
		return nil, fmt.Errorf("linear transformation has invalid N1=%d", matrix.N1)
	}
	colsLog := matrix.LogDimensions.Cols
	if colsLog < 0 || colsLog >= params.LogN() {
		return nil, fmt.Errorf("linear transformation has invalid LogDimensions.Cols=%d", colsLog)
	}
	cols := 1 << uint(colsLog)
	mirrored := &FastStorageLinearTransformation{
		params: params, width: storageWidth, levelQ: matrix.LevelQ, scale: cloneFastScale(matrix.Scale),
		metadata: cloneFastPlaintextMetadata(matrix.PlaintextMetaData), n1: matrix.N1,
		cols: cols, diagonals: make(map[int]*FastStoragePlaintext, len(matrix.Vec)),
		sumL1Norm: new(big.Int), maxAbs: new(big.Int), maxL1Norm: new(big.Int),
	}
	for diagonal, qp := range matrix.Vec {
		key := diagonal & (cols - 1)
		if _, exists := mirrored.diagonals[key]; exists {
			return nil, fmt.Errorf("linear transformation has duplicate diagonal index modulo %d: %d", cols, diagonal)
		}
		pt, err := NewFastStoragePlaintextMirrorWithWidth(params, qp.Q, matrix.LevelQ, matrix.PlaintextMetaData, matrix.IsNTT, matrix.IsMontgomery, storageWidth)
		if err != nil {
			return nil, fmt.Errorf("mirror diagonal %d: %w", diagonal, err)
		}
		mirrored.diagonals[key] = pt
		mirrored.sumL1Norm.Add(mirrored.sumL1Norm, pt.l1Norm)
		if pt.maxAbs.Cmp(mirrored.maxAbs) > 0 {
			mirrored.maxAbs.Set(pt.maxAbs)
		}
		if pt.l1Norm.Cmp(mirrored.maxL1Norm) > 0 {
			mirrored.maxL1Norm.Set(pt.l1Norm)
		}
	}
	if matrix.N1 != 0 {
		mirrored.index, _, _ = matrix.BSGSIndex()
		if len(mirrored.index) == 0 {
			return nil, errors.New("linear transformation declares BSGS but has an empty BSGS schedule")
		}
	}
	return mirrored, nil
}

func (matrix *FastStorageLinearTransformation) DiagonalCount() int {
	if matrix == nil {
		return 0
	}
	return len(matrix.diagonals)
}

func (matrix *FastStorageLinearTransformation) LevelQ() int {
	if matrix == nil {
		return -1
	}
	return matrix.levelQ
}

func (matrix *FastStorageLinearTransformation) Scale() rlwe.Scale {
	if matrix == nil {
		return rlwe.Scale{}
	}
	return cloneFastScale(matrix.scale)
}

func (matrix *FastStorageLinearTransformation) N1() int {
	if matrix == nil {
		return 0
	}
	return matrix.n1
}

func (matrix *FastStorageLinearTransformation) MaxAbsCoefficient() *big.Int {
	if matrix == nil || matrix.maxAbs == nil {
		return nil
	}
	return new(big.Int).Set(matrix.maxAbs)
}

func (matrix *FastStorageLinearTransformation) MaxDiagonalL1Norm() *big.Int {
	if matrix == nil || matrix.maxL1Norm == nil {
		return nil
	}
	return new(big.Int).Set(matrix.maxL1Norm)
}

func (matrix *FastStorageLinearTransformation) SumDiagonalL1Norm() *big.Int {
	if matrix == nil || matrix.sumL1Norm == nil {
		return nil
	}
	return new(big.Int).Set(matrix.sumL1Norm)
}

// FastStorageLinearTransform evaluates a mirrored diagonal transform in
// private F. It uses the original matrix's direct or BSGS schedule and
// preflights the conservative L1 capacity bound before arithmetic begins.
func FastStorageLinearTransform(ct *FastCiphertext, matrix *FastStorageLinearTransformation) (*FastCiphertext, error) {
	if err := validateFastStorageState(ct); err != nil {
		return nil, fmt.Errorf("ciphertext: %w", err)
	}
	if matrix == nil || len(matrix.diagonals) == 0 {
		return nil, errors.New("private-F linear transformation cannot be nil or empty")
	}
	if ct.activeStorageWidth != matrix.width || ct.Degree() != 1 || !ct.metadata.IsNTT || ct.metadata.IsMontgomery {
		return nil, fmt.Errorf("private-F LinearTransform requires matching-width degree-one NTT ordinary-residue storage (ciphertext=%d matrix=%d)", ct.activeStorageWidth, matrix.width)
	}
	if ct.N() != matrix.params.N() || !ct.params.Equal(&matrix.params) || ct.logicalLevel > matrix.levelQ {
		return nil, errors.New("private-F LinearTransform CKKS parameters or LevelQ are incompatible with input")
	}
	if ct.metadata.IsBatched != matrix.metadata.IsBatched {
		return nil, errors.New("private-F LinearTransform batching metadata is incompatible with input")
	}
	bounds := make([]*big.Int, len(ct.componentBounds))
	for component, bound := range ct.componentBounds {
		bounds[component] = new(big.Int).Mul(bound, matrix.sumL1Norm)
	}
	if err := validateBoundsFit(ct.basis, bounds, ct.activeStorageWidth); err != nil {
		return nil, fmt.Errorf("private-F LinearTransform capacity preflight: %w", err)
	}

	output := newFastCiphertextWithBasis(ct.params, 1, ct.logicalLevel, ct.activeStorageWidth, ct.basis)
	output.metadata = cloneFastMetadata(ct.metadata)
	output.metadata.Scale = ct.metadata.Scale.Mul(matrix.scale)
	output.componentBounds = bounds
	if matrix.n1 == 0 {
		if err := fastStorageLinearTransformDirect(ct, matrix, output); err != nil {
			return nil, err
		}
	} else if err := fastStorageLinearTransformBSGS(ct, matrix, output); err != nil {
		return nil, err
	}
	if err := validateFastStorageState(output); err != nil {
		return nil, fmt.Errorf("invalid private-F LinearTransform output: %w", err)
	}
	return output, nil
}

func fastStorageLinearTransformDirect(input *FastCiphertext, matrix *FastStorageLinearTransformation, output *FastCiphertext) error {
	keys := make([]int, 0, len(matrix.diagonals))
	for key := range matrix.diagonals {
		keys = append(keys, key)
	}
	sort.Ints(keys)
	for _, diagonal := range keys {
		pt := matrix.diagonals[diagonal]
		rotated, err := fastStorageRotatedComponents(input, diagonal)
		if err != nil {
			return fmt.Errorf("private-F diagonal %d automorphism: %w", diagonal, err)
		}
		if err := fastStorageMulAddTerm(output, rotated, pt); err != nil {
			return fmt.Errorf("private-F diagonal %d multiply: %w", diagonal, err)
		}
	}
	return nil
}

func fastStorageLinearTransformBSGS(input *FastCiphertext, matrix *FastStorageLinearTransformation, output *FastCiphertext) error {
	babies := make(map[int][]ring.Poly)
	babyKeys := make(map[int]bool)
	for _, babyIndexes := range matrix.index {
		for _, i := range babyIndexes {
			babyKeys[i] = true
		}
	}
	orderedBabies := make([]int, 0, len(babyKeys))
	for i := range babyKeys {
		orderedBabies = append(orderedBabies, i)
	}
	sort.Ints(orderedBabies)
	for _, i := range orderedBabies {
		rotated, err := fastStorageRotatedComponents(input, i)
		if err != nil {
			return fmt.Errorf("private-F baby automorphism %d: %w", i, err)
		}
		babies[i] = rotated
	}

	outerKeys := make([]int, 0, len(matrix.index))
	for j := range matrix.index {
		outerKeys = append(outerKeys, j)
	}
	sort.Ints(outerKeys)
	for _, j := range outerKeys {
		inner := newFastStorageComponentPair(input.N(), input.activeStorageWidth)
		babyIndexes := append([]int(nil), matrix.index[j]...)
		sort.Ints(babyIndexes)
		for _, i := range babyIndexes {
			diagonal := (j + i) % matrix.cols
			pt, ok := matrix.diagonals[diagonal]
			if !ok {
				return fmt.Errorf("private-F BSGS schedule references missing diagonal %d", diagonal)
			}
			if err := fastStorageMulAddTermPair(inner, babies[i], pt, input.basis); err != nil {
				return fmt.Errorf("private-F BSGS diagonal %d multiply: %w", diagonal, err)
			}
		}
		rotated, err := fastStorageRotatePair(inner, input, j)
		if err != nil {
			return fmt.Errorf("private-F giant automorphism %d: %w", j, err)
		}
		addFastStoragePair(output, rotated, input.basis)
	}
	return nil
}

func newFastStorageComponentPair(N, width int) []ring.Poly {
	return []ring.Poly{ring.NewPoly(N, width-1), ring.NewPoly(N, width-1)}
}

func fastStorageRotatedComponents(input *FastCiphertext, rotation int) ([]ring.Poly, error) {
	rotated := newFastStorageComponentPair(input.N(), input.activeStorageWidth)
	if rotation == 0 {
		for component := 0; component < 2; component++ {
			for row := 0; row < input.activeStorageWidth; row++ {
				copy(rotated[component].Coeffs[row], input.value[component].Coeffs[row])
			}
		}
		return rotated, nil
	}
	index, err := ring.AutomorphismNTTIndex(input.N(), input.params.RingQ().NthRoot(), input.params.GaloisElement(rotation))
	if err != nil {
		return nil, err
	}
	for component := 0; component < 2; component++ {
		for row := 0; row < input.activeStorageWidth; row++ {
			permuteStorageNTT(input.value[component].Coeffs[row], rotated[component].Coeffs[row], index)
		}
	}
	return rotated, nil
}

func fastStorageRotatePair(input []ring.Poly, ct *FastCiphertext, rotation int) ([]ring.Poly, error) {
	if rotation == 0 {
		return input, nil
	}
	index, err := ring.AutomorphismNTTIndex(ct.N(), ct.params.RingQ().NthRoot(), ct.params.GaloisElement(rotation))
	if err != nil {
		return nil, err
	}
	output := newFastStorageComponentPair(ct.N(), ct.activeStorageWidth)
	for component := 0; component < 2; component++ {
		for row := 0; row < ct.activeStorageWidth; row++ {
			permuteStorageNTT(input[component].Coeffs[row], output[component].Coeffs[row], index)
		}
	}
	return output, nil
}

func fastStorageMulAddTerm(output *FastCiphertext, rotated []ring.Poly, pt *FastStoragePlaintext) error {
	for row := 0; row < output.activeStorageWidth; row++ {
		subring, _ := output.basis.subring(row)
		for component := 0; component < 2; component++ {
			subring.MulCoeffsBarrettThenAdd(rotated[component].Coeffs[row], pt.value.Coeffs[row], output.value[component].Coeffs[row])
		}
	}
	return nil
}

func fastStorageMulAddTermPair(output, rotated []ring.Poly, pt *FastStoragePlaintext, basis fastStorageBasis) error {
	for row := 0; row < len(output[0].Coeffs); row++ {
		subring, _ := basis.subring(row)
		for component := 0; component < 2; component++ {
			subring.MulCoeffsBarrettThenAdd(rotated[component].Coeffs[row], pt.value.Coeffs[row], output[component].Coeffs[row])
		}
	}
	return nil
}

func addFastStoragePair(output *FastCiphertext, input []ring.Poly, basis fastStorageBasis) {
	for row := 0; row < output.activeStorageWidth; row++ {
		subring, _ := basis.subring(row)
		for component := 0; component < 2; component++ {
			subring.Add(output.value[component].Coeffs[row], input[component].Coeffs[row], output.value[component].Coeffs[row])
		}
	}
}
