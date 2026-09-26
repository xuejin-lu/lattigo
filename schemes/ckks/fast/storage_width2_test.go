package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	commonlintrans "github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func TestFastStorageContractThreeToTwoExactnessAndMetadata(t *testing.T) {
	params := storageBoundaryParameters(t)
	polynomials := [][]int64{{0, 1, -1, 17, -23}, {8, -5, 2}}
	for _, ntt := range []bool{false, true} {
		name := "coeff"
		if ntt {
			name = "ntt"
		}
		t.Run(name, func(t *testing.T) {
			input := newSyntheticStorageCiphertext(t, params, polynomials, 3, 3, ntt, nil)
			input.metadata.Scale = rlwe.NewScale(123)
			input.metadata.LogDimensions = ring.Dimensions{Cols: 2}
			before := snapshotFastStorageCiphertext(input)
			want := decodeFastStoragePolynomials(t, input)

			output, err := input.ContractStorage(2)
			require.NoError(t, err)
			require.Equal(t, 2, output.StorageWidth())
			require.Equal(t, input.LogicalLevel(), output.LogicalLevel())
			require.Equal(t, input.Degree(), output.Degree())
			require.Equal(t, input.MetaData(), output.MetaData())
			require.Equal(t, input.ComponentBounds(), output.ComponentBounds())
			requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
			require.NoError(t, validateFastStorageState(output))
			require.Equal(t, before, snapshotFastStorageCiphertext(input), "contraction must not mutate its input")
		})
	}
}

func TestFastStorageContractCapacityFailureIsTransactional(t *testing.T) {
	params := storageBoundaryParameters(t)
	basis, err := fastStorageBasisForLogN(params.LogN())
	require.NoError(t, err)
	capacity2, err := basis.centeredCapacity(2)
	require.NoError(t, err)
	tooLarge := new(big.Int).Add(capacity2, big.NewInt(1))
	input := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{tooLarge}}, 3, 3, true, []*big.Int{tooLarge})
	before := snapshotFastStorageCiphertext(input)
	_, err = input.ContractStorage(2)
	require.ErrorContains(t, err, "storage contraction capacity preflight")
	require.Equal(t, before, snapshotFastStorageCiphertext(input))
	require.NoError(t, validateFastStorageState(input))
}

func TestFastStorageWidth2PlaintextMirrorAndArithmeticExactness(t *testing.T) {
	params := fastStoragePlaintextTestParameters(t)
	level := 2
	values := make([]*big.Int, params.N())
	for i := range values {
		values[i] = new(big.Int)
	}
	values[0].SetInt64(2)
	values[1].SetInt64(-3)
	values[4].SetInt64(5)
	metadata := rlwe.PlaintextMetaData{Scale: rlwe.NewScale(32), IsBatched: true, LogDimensions: ring.Dimensions{Cols: 2}}
	logical := logicalQPolynomialFromIntegers(params, level, values)
	pt2, err := NewFastStoragePlaintextMirrorWithWidth(params, logical, level, metadata, true, true, 2)
	require.NoError(t, err)
	pt3, err := NewFastStoragePlaintextMirrorWithWidth(params, logical, level, metadata, true, true, 3)
	require.NoError(t, err)
	require.Equal(t, 2, pt2.StorageWidth())
	require.Equal(t, 3, pt3.StorageWidth())
	requireBigIntVectorsEqual(t, decodeFastStoragePlaintext(t, pt3), decodeFastStoragePlaintext(t, pt2))
	require.Equal(t, pt3.MaxAbsCoefficient(), pt2.MaxAbsCoefficient())
	require.Equal(t, pt3.L1Norm(), pt2.L1Norm())
	require.True(t, pt3.Scale().Equal(pt2.Scale()))
	basis, err := fastStorageBasisForLogN(params.LogN())
	require.NoError(t, err)
	capacity2, err := basis.centeredCapacity(2)
	require.NoError(t, err)
	tooLarge := make([]*big.Int, params.N())
	for i := range tooLarge {
		tooLarge[i] = new(big.Int)
	}
	tooLarge[0].Add(capacity2, big.NewInt(1))
	tooLargeLogical := logicalQPolynomialFromIntegers(params, level, tooLarge)
	_, err = NewFastStoragePlaintextMirrorWithWidth(params, tooLargeLogical, level, metadata, true, true, 3)
	require.NoError(t, err, "the coefficient remains within width-3 capacity")
	_, err = NewFastStoragePlaintextMirrorWithWidth(params, tooLargeLogical, level, metadata, true, true, 2)
	require.ErrorContains(t, err, "strict width-2 centered capacity")

	inputPolys := [][]*big.Int{{big.NewInt(3), big.NewInt(-2), big.NewInt(5)}, {big.NewInt(-7), big.NewInt(4)}}
	ct3 := newFastStorageCiphertextFromBig(t, params, inputPolys, level, 3, true, nil)
	ct3.metadata.LogDimensions = metadata.LogDimensions
	ct2, err := ct3.ContractStorage(2)
	require.NoError(t, err)

	mul2, err := FastStorageMulPlaintext(ct2, pt2)
	require.NoError(t, err)
	mul3, err := FastStorageMulPlaintext(ct3, pt3)
	require.NoError(t, err)
	plainLift := decodeFastStoragePlaintext(t, pt2)
	wantMul := make([][]*big.Int, len(inputPolys))
	for component := range inputPolys {
		wantMul[component] = independentNegacyclicProduct(inputPolys[component], plainLift, params.N())
	}
	requireBigIntMatricesEqual(t, wantMul, decodeFastStoragePolynomials(t, mul2))
	requireBigIntMatricesEqual(t, decodeFastStoragePolynomials(t, mul3), decodeFastStoragePolynomials(t, mul2))
	require.Equal(t, mul3.ComponentBounds(), mul2.ComponentBounds())

	galEl := params.GaloisElement(3)
	auto2, err := FastStorageAutomorphism(ct2, galEl)
	require.NoError(t, err)
	auto3, err := FastStorageAutomorphism(ct3, galEl)
	require.NoError(t, err)
	wantAuto := make([][]*big.Int, len(inputPolys))
	for component := range inputPolys {
		wantAuto[component] = independentAutomorphism(inputPolys[component], params.N(), galEl)
	}
	requireBigIntMatricesEqual(t, wantAuto, decodeFastStoragePolynomials(t, auto2))
	requireBigIntMatricesEqual(t, decodeFastStoragePolynomials(t, auto3), decodeFastStoragePolynomials(t, auto2))

	diagonals := map[int][]complex128{
		0: {2, -1, 0, 3},
		1: {1, 0, 2, -2},
		3: {-1, 2, 1, 0},
	}
	for _, ratio := range []int{-1, 1} {
		logicalMatrix := encodedTestLinearTransformation(t, params, level, diagonals, ratio)
		matrix2, err := NewFastStorageLinearTransformationWithWidth(params, commonlintrans.LinearTransformation(logicalMatrix), 2)
		require.NoError(t, err)
		matrix3, err := NewFastStorageLinearTransformationWithWidth(params, commonlintrans.LinearTransformation(logicalMatrix), 3)
		require.NoError(t, err)
		out2, err := FastStorageLinearTransform(ct2, matrix2)
		require.NoError(t, err)
		out3, err := FastStorageLinearTransform(ct3, matrix3)
		require.NoError(t, err)
		requireBigIntMatricesEqual(t, decodeFastStoragePolynomials(t, out3), decodeFastStoragePolynomials(t, out2))
		require.Equal(t, out3.ComponentBounds(), out2.ComponentBounds())
	}
}

func requireBigIntVectorsEqual(t *testing.T, expected, actual []*big.Int) {
	t.Helper()
	require.Len(t, actual, len(expected))
	for i := range expected {
		require.Zero(t, expected[i].Cmp(actual[i]), "coefficient=%d expected=%s actual=%s", i, expected[i], actual[i])
	}
}
