package fast

import (
	"math/big"
	"strconv"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFastStorageModUpLevel0CanonicalCenteredQ0(t *testing.T) {
	params := storageBoundaryParameters(t)
	q0 := new(big.Int).SetUint64(params.Q()[0])
	half := new(big.Int).Rsh(new(big.Int).Set(q0), 1)
	values := []*big.Int{
		big.NewInt(0), big.NewInt(1), big.NewInt(-1),
		new(big.Int).Set(q0), new(big.Int).Neg(new(big.Int).Set(q0)),
		new(big.Int).Lsh(new(big.Int).Set(q0), 1), new(big.Int).Neg(new(big.Int).Lsh(new(big.Int).Set(q0), 1)),
		new(big.Int).Set(half), new(big.Int).Add(new(big.Int).Set(half), big.NewInt(1)),
		new(big.Int).Neg(new(big.Int).Set(half)), new(big.Int).Neg(new(big.Int).Add(new(big.Int).Set(half), big.NewInt(1))),
		new(big.Int).Add(big.NewInt(3), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(2))),
		new(big.Int).Add(big.NewInt(-7), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(-3))),
	}
	input := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{values}, 0, 3, false, nil)
	output, err := FastStorageModUpLevel0(input, 1)
	require.NoError(t, err)
	want := make([][]*big.Int, 1)
	want[0] = make([]*big.Int, params.N())
	for i, value := range values {
		want[0][i] = centeredBigIntModulo(value, q0)
	}
	for i := len(values); i < params.N(); i++ {
		want[0][i] = new(big.Int)
	}
	requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
}

func TestFastStorageModUpLevel0CoefficientExactnessAndBounds(t *testing.T) {
	params := storageBoundaryParameters(t)
	q0 := new(big.Int).SetUint64(params.Q()[0])
	values := [][]*big.Int{
		{
			new(big.Int).Add(new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(113)), big.NewInt(9)),
			new(big.Int).Add(new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(-71)), big.NewInt(-23)),
			new(big.Int).Add(new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(4)), new(big.Int).Rsh(new(big.Int).Set(q0), 1)),
		},
		{
			new(big.Int).Add(new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(-8)), big.NewInt(31)),
			new(big.Int).Add(new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(19)), big.NewInt(-47)),
		},
		{new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(33)), big.NewInt(-1)},
	}
	input := newFastStorageCiphertextFromBig(t, params, values, 0, 3, false, nil)
	output, err := FastStorageModUpLevel0(input, 2)
	require.NoError(t, err)
	want := padBigIntMatrixToN(centeredMatrixModuloQ0(values, q0), params.N())
	requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
	require.Equal(t, 3, output.StorageWidth())
	for component := range values {
		wantBound := maxBigIntMagnitude(want[component])
		require.Zero(t, wantBound.Cmp(output.ComponentBounds()[component]), "component=%d", component)
		require.Less(t, output.ComponentBounds()[component].Cmp(input.ComponentBounds()[component]), 0,
			"canonicalization must discard the large q0-multiple contribution, component=%d", component)
	}
}

func TestFastStorageModUpLevel0NTTExactnessAndDomainPreservation(t *testing.T) {
	params := storageBoundaryParameters(t)
	q0 := new(big.Int).SetUint64(params.Q()[0])
	values := [][]*big.Int{
		{new(big.Int).Add(big.NewInt(5), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(7))), big.NewInt(-13)},
		{new(big.Int).Sub(big.NewInt(-2), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(11))), new(big.Int).Set(q0)},
	}
	want := padBigIntMatrixToN(centeredMatrixModuloQ0(values, q0), params.N())
	for _, ntt := range []bool{false, true} {
		input := newFastStorageCiphertextFromBig(t, params, values, 0, 3, ntt, nil)
		output, err := FastStorageModUpLevel0(input, 3)
		require.NoError(t, err)
		require.Equal(t, ntt, output.IsNTT())
		require.False(t, output.IsMontgomery())
		requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
	}
}

func TestFastStorageModUpLevel0TargetLevelsAndLogicalExportOracle(t *testing.T) {
	params := storageBoundaryParameters(t)
	q0 := new(big.Int).SetUint64(params.Q()[0])
	values := [][]*big.Int{
		{new(big.Int).Add(big.NewInt(3), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(2))), new(big.Int).Sub(big.NewInt(-4), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(3)))},
		{new(big.Int).Add(big.NewInt(17), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(-5))), new(big.Int).Neg(new(big.Int).Set(q0))},
	}
	for _, targetLevel := range []int{1, 3, params.MaxLevel()} {
		for _, ntt := range []bool{false, true} {
			t.Run(fmtModUpCase(targetLevel, ntt), func(t *testing.T) {
				input := newFastStorageCiphertextFromBig(t, params, values, 0, 3, ntt, nil)
				output, err := FastStorageModUpLevel0(input, targetLevel)
				require.NoError(t, err)
				require.Equal(t, targetLevel, output.LogicalLevel())
				require.Equal(t, 3, output.StorageWidth())
				require.Equal(t, ntt, output.IsNTT())
				assertFastStorageModUpLogicalOracle(t, input, output)
				for component := range output.value {
					require.Len(t, output.value[component].Coeffs, 3, "private-F container must retain exactly three rows")
				}
				logical, err := output.ExportToLogical(FastCiphertextDomain{})
				require.NoError(t, err)
				require.Equal(t, targetLevel+1, len(logical.Value[0].Coeffs), "logical export must materialize exactly the requested q rows")
			})
		}
	}
}

func TestFastStorageModUpLevel0PreservesScaleMetadataDegreeAndParameters(t *testing.T) {
	params := storageBoundaryParameters(t)
	paramsBefore := params.ParametersLiteral()
	values := [][]*big.Int{{big.NewInt(3), big.NewInt(-7)}, {big.NewInt(11)}, {big.NewInt(-13)}}
	input := newFastStorageCiphertextFromBig(t, params, values, 0, 3, true, nil)
	input.metadata.Scale = rlwe.NewScale(123456789012345)
	input.metadata.IsBatched = false
	input.metadata.IsBitReversed = true
	input.metadata.LogDimensions = ring.Dimensions{Rows: 2, Cols: 3}
	beforeMetadata := input.MetaData()
	inputRows := cloneStorageRows(input.value)
	inputBounds := input.ComponentBounds()

	output, err := FastStorageModUpLevel0(input, params.MaxLevel())
	require.NoError(t, err)
	require.True(t, beforeMetadata.Scale.Equal(output.Scale()))
	require.Equal(t, input.Degree(), output.Degree())
	require.Equal(t, 3, output.StorageWidth())
	require.Equal(t, input.IsNTT(), output.IsNTT())
	require.False(t, output.IsMontgomery())
	gotMetadata := output.MetaData()
	require.Equal(t, beforeMetadata.IsBatched, gotMetadata.IsBatched)
	require.Equal(t, beforeMetadata.IsBitReversed, gotMetadata.IsBitReversed)
	require.Equal(t, beforeMetadata.LogDimensions, gotMetadata.LogDimensions)
	outputParameters := output.Parameters()
	require.True(t, params.Equal(&outputParameters))
	require.Equal(t, paramsBefore, outputParameters.ParametersLiteral())
	require.Equal(t, paramsBefore, params.ParametersLiteral(), "ModUp must not mutate public parameters")
	require.Equal(t, inputRows, input.value, "input rows must remain unchanged")
	require.Equal(t, inputBounds, input.ComponentBounds(), "input bounds must remain unchanged")
}

func TestFastStorageModUpLevel0NoncanonicalLiftTheorem(t *testing.T) {
	params := storageBoundaryParameters(t)
	q0 := new(big.Int).SetUint64(params.Q()[0])
	canonical := [][]*big.Int{{big.NewInt(3), big.NewInt(-4), big.NewInt(7)}, {big.NewInt(-2), big.NewInt(5)}}
	makeLift := func(ks [][]int64) [][]*big.Int {
		values := make([][]*big.Int, len(canonical))
		for component := range canonical {
			values[component] = make([]*big.Int, len(canonical[component]))
			for coefficient, value := range canonical[component] {
				values[component][coefficient] = new(big.Int).Add(new(big.Int).Set(value), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(ks[component][coefficient])))
			}
		}
		return values
	}
	liftA := makeLift([][]int64{{2, -3, 5}, {-7, 11}})
	liftB := makeLift([][]int64{{-13, 17, -19}, {23, -29}})
	for _, ntt := range []bool{false, true} {
		inputA := newFastStorageCiphertextFromBig(t, params, liftA, 0, 3, ntt, nil)
		inputB := newFastStorageCiphertextFromBig(t, params, liftB, 0, 3, ntt, nil)
		outputA, err := FastStorageModUpLevel0(inputA, 3)
		require.NoError(t, err)
		outputB, err := FastStorageModUpLevel0(inputB, 3)
		require.NoError(t, err)
		require.Equal(t, outputA.value, outputB.value, "different c+kq0 lifts must canonicalize identically")
		require.Equal(t, padBigIntMatrixToN(centeredMatrixModuloQ0(liftA, q0), params.N()), decodeFastStoragePolynomials(t, outputA))
		logicalA, err := outputA.ExportToLogical(FastCiphertextDomain{})
		require.NoError(t, err)
		logicalB, err := outputB.ExportToLogical(FastCiphertextDomain{})
		require.NoError(t, err)
		require.Equal(t, logicalA.Value, logicalB.Value)
	}
}

func TestFastStorageModUpLevel0FailureIsTransactional(t *testing.T) {
	params := storageBoundaryParameters(t)
	tests := []struct {
		name        string
		width       int
		targetLevel int
		mutate      func(*FastCiphertext)
		wantError   string
	}{
		{name: "nonzero-level", width: 3, targetLevel: 1, mutate: func(ct *FastCiphertext) { ct.logicalLevel = 1 }, wantError: "requires logical Level 0"},
		{name: "width-one", width: 1, targetLevel: 1, wantError: "requires storage width 3"},
		{name: "width-two", width: 2, targetLevel: 1, wantError: "requires storage width 3"},
		{name: "target-zero", width: 3, targetLevel: 0, wantError: "outside [1,"},
		{name: "target-negative", width: 3, targetLevel: -1, wantError: "outside [1,"},
		{name: "target-too-high", width: 3, targetLevel: params.MaxLevel() + 1, wantError: "outside [1,"},
		{name: "montgomery", width: 3, targetLevel: 1, mutate: func(ct *FastCiphertext) { ct.metadata.IsMontgomery = true }, wantError: "Montgomery"},
		{name: "malformed-row", width: 3, targetLevel: 1, mutate: func(ct *FastCiphertext) { ct.value[0].Coeffs[1] = ct.value[0].Coeffs[1][:1] }, wantError: "expected N"},
		{name: "malformed-bound", width: 3, targetLevel: 1, mutate: func(ct *FastCiphertext) { ct.componentBounds[0] = nil }, wantError: "bound cannot be nil"},
		{name: "non-positive-scale", width: 3, targetLevel: 1, mutate: func(ct *FastCiphertext) { ct.metadata.Scale = rlwe.NewScale(0) }, wantError: "scale must be positive"},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			input := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{big.NewInt(7)}}, 0, test.width, false, nil)
			if test.mutate != nil {
				test.mutate(input)
			}
			before := snapshotFastStorageCiphertext(input)
			output, err := FastStorageModUpLevel0(input, test.targetLevel)
			require.Nil(t, output)
			require.ErrorContains(t, err, test.wantError)
			require.Equal(t, before, snapshotFastStorageCiphertext(input), "failed ModUp must not mutate the input")
		})
	}
	output, err := FastStorageModUpLevel0(nil, 1)
	require.Nil(t, output)
	require.ErrorContains(t, err, "cannot be nil")

	ciParams, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 4, LogQ: []int{50, 50}, RingType: ring.ConjugateInvariant})
	require.NoError(t, err)
	_, err = NewFastCiphertext(ciParams, 1, 0, 3)
	require.ErrorContains(t, err, "Standard ring")
}

func TestFastStorageModUpLevel0RescaleBoundaryChain(t *testing.T) {
	params := storageBoundaryParameters(t)
	q0 := new(big.Int).SetUint64(params.Q()[0])
	q1 := new(big.Int).SetUint64(params.Q()[1])
	large := new(big.Int).Add(q0, big.NewInt(17))
	inputValue := new(big.Int).Mul(new(big.Int).Set(q1), large)
	input := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{inputValue}, {new(big.Int).Neg(new(big.Int).Set(inputValue))}}, 1, 3, true, nil)
	input.metadata.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 100))

	rescaled, err := FastStorageRescale(input)
	require.NoError(t, err)
	require.Equal(t, 0, rescaled.LogicalLevel())
	rescaledValues := decodeFastStoragePolynomials(t, rescaled)
	require.Equal(t, q0.Int64()+17, rescaledValues[0][0].Int64(), "fixture must produce a noncanonical q0 lift")
	rescaledScale := rescaled.Scale()
	output, err := FastStorageModUpLevel0(rescaled, params.MaxLevel())
	require.NoError(t, err)
	require.Equal(t, params.MaxLevel(), output.LogicalLevel())
	require.True(t, rescaledScale.Equal(output.Scale()), "ModUp must leave the Rescale scale unchanged")
	want := centeredMatrixModuloQ0(rescaledValues, q0)
	requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
	assertFastStorageModUpLogicalOracle(t, rescaled, output)
}

func assertFastStorageModUpLogicalOracle(t *testing.T, input, output *FastCiphertext) {
	t.Helper()
	inputLogical, err := input.ExportToLogical(FastCiphertextDomain{})
	require.NoError(t, err)
	outputLogical, err := output.ExportToLogical(FastCiphertextDomain{})
	require.NoError(t, err)
	q0 := input.params.Q()[0]
	for component := range inputLogical.Value {
		for coefficient, residue := range inputLogical.Value[component].Coeffs[0] {
			canonical := centeredStorageResidue(residue, q0)
			value := signedStorageBig(canonical)
			for logicalIndex := 0; logicalIndex <= output.LogicalLevel(); logicalIndex++ {
				modulus := new(big.Int).SetUint64(input.params.Q()[logicalIndex])
				want := new(big.Int).Mod(new(big.Int).Set(value), modulus).Uint64()
				require.Equal(t, want, outputLogical.Value[component].Coeffs[logicalIndex][coefficient],
					"component=%d coefficient=%d logical-index=%d", component, coefficient, logicalIndex)
			}
		}
	}
}

func centeredBigIntModulo(value, modulus *big.Int) *big.Int {
	residue := new(big.Int).Mod(new(big.Int).Set(value), modulus)
	if residue.Cmp(new(big.Int).Rsh(new(big.Int).Set(modulus), 1)) > 0 {
		residue.Sub(residue, modulus)
	}
	return residue
}

func centeredMatrixModuloQ0(values [][]*big.Int, q0 *big.Int) [][]*big.Int {
	want := make([][]*big.Int, len(values))
	for component := range values {
		want[component] = make([]*big.Int, len(values[component]))
		for coefficient, value := range values[component] {
			want[component][coefficient] = centeredBigIntModulo(value, q0)
		}
	}
	return want
}

func padBigIntMatrixToN(values [][]*big.Int, n int) [][]*big.Int {
	padded := make([][]*big.Int, len(values))
	for component := range values {
		padded[component] = make([]*big.Int, n)
		for coefficient := range padded[component] {
			padded[component][coefficient] = new(big.Int)
			if coefficient < len(values[component]) {
				padded[component][coefficient].Set(values[component][coefficient])
			}
		}
	}
	return padded
}

func maxBigIntMagnitude(values []*big.Int) *big.Int {
	maximum := new(big.Int)
	for _, value := range values {
		magnitude := new(big.Int).Abs(new(big.Int).Set(value))
		if magnitude.Cmp(maximum) > 0 {
			maximum.Set(magnitude)
		}
	}
	return maximum
}

func fmtModUpCase(level int, ntt bool) string {
	domain := "coeff"
	if ntt {
		domain = "ntt"
	}
	return "level-" + strconv.Itoa(level) + "/" + domain
}
