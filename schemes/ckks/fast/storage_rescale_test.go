package fast

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestRoundStorageIntegerByUint64SignedBoundaries(t *testing.T) {
	params := storageBoundaryParameters(t)
	q := params.Q()[3]
	qBig := new(big.Int).SetUint64(q)
	values := []*big.Int{
		big.NewInt(0), big.NewInt(1), big.NewInt(-1),
		new(big.Int).Set(qBig), new(big.Int).Neg(new(big.Int).Set(qBig)),
		new(big.Int).Rsh(new(big.Int).Sub(new(big.Int).Set(qBig), big.NewInt(1)), 1),
		new(big.Int).Neg(new(big.Int).Rsh(new(big.Int).Sub(new(big.Int).Set(qBig), big.NewInt(1)), 1)),
		new(big.Int).Rsh(new(big.Int).Add(new(big.Int).Set(qBig), big.NewInt(1)), 1),
		new(big.Int).Neg(new(big.Int).Rsh(new(big.Int).Add(new(big.Int).Set(qBig), big.NewInt(1)), 1)),
	}
	for k := int64(1); k <= 7; k++ {
		for _, delta := range []int64{-2, -1, 0, 1, 2} {
			value := new(big.Int).Mul(qBig, big.NewInt(k))
			value.Add(value, big.NewInt(delta))
			values = append(values, value, new(big.Int).Neg(new(big.Int).Set(value)))
		}
	}
	for _, value := range values {
		got, err := roundStorageIntegerByUint64(storageIntegerFromBig(value), q)
		require.NoError(t, err)
		require.Zero(t, roundBigIntByUint64(value, q).Cmp(signedStorageBig(got)), "value=%s q=%d", value, q)
	}
	_, err := roundStorageIntegerByUint64(storageInteger{}, 0)
	require.ErrorContains(t, err, "cannot be zero")
	_, err = roundStorageIntegerByUint64(storageInteger{}, 2)
	require.ErrorContains(t, err, "odd divisor")
}

func TestRoundStorageIntegerByUint64AgainstBigIntAcross192Bits(t *testing.T) {
	paramsA := storageBoundaryParameters(t)
	paramsB := storageBoundaryParametersLogN(t, 5)
	divisors := append(paramsA.Q()[1:4], paramsB.Q()[1:3]...)
	require.NotEmpty(t, divisors)
	rng := rand.New(rand.NewSource(4004))
	for width := 1; width <= 3; width++ {
		basis, err := newFastStorageBasis(4)
		require.NoError(t, err)
		capacity, err := basis.centeredCapacity(width)
		require.NoError(t, err)
		magnitudes := []*big.Int{new(big.Int), big.NewInt(1), new(big.Int).Set(capacity)}
		if capacity.Sign() > 0 {
			magnitudes = append(magnitudes, new(big.Int).Sub(new(big.Int).Set(capacity), big.NewInt(1)))
		}
		for i := 0; i < 40; i++ {
			value := new(big.Int).SetUint64(rng.Uint64())
			value.Lsh(value, 64).Or(value, new(big.Int).SetUint64(rng.Uint64()))
			value.Lsh(value, 64).Or(value, new(big.Int).SetUint64(rng.Uint64()))
			value.Mod(value, new(big.Int).Add(new(big.Int).Set(capacity), big.NewInt(1)))
			magnitudes = append(magnitudes, value)
		}
		for _, divisor := range divisors {
			for index, magnitude := range magnitudes {
				value := new(big.Int).Set(magnitude)
				if index&1 == 1 {
					value.Neg(value)
				}
				got, err := roundStorageIntegerByUint64(storageIntegerFromBig(value), divisor)
				require.NoError(t, err)
				require.Zero(t, roundBigIntByUint64(value, divisor).Cmp(signedStorageBig(got)), "width=%d divisor=%d value=%s", width, divisor, value)
			}
		}
	}
	max192 := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), 192), big.NewInt(1))
	for _, divisor := range divisors {
		for _, value := range []*big.Int{max192, new(big.Int).Neg(new(big.Int).Set(max192))} {
			got, err := roundStorageIntegerByUint64(storageIntegerFromBig(value), divisor)
			require.NoError(t, err)
			require.Zero(t, roundBigIntByUint64(value, divisor).Cmp(signedStorageBig(got)))
		}
	}
}

func TestFastStorageRescaleUsesLogicalDivisorAcrossStorageWidths(t *testing.T) {
	params := storageBoundaryParameters(t)
	level := 3
	q := params.Q()[level]
	values := [][]*big.Int{{
		new(big.Int).Add(new(big.Int).Mul(new(big.Int).SetUint64(q), big.NewInt(3)), big.NewInt(2)),
		new(big.Int).Neg(new(big.Int).Add(new(big.Int).Mul(new(big.Int).SetUint64(q), big.NewInt(2)), big.NewInt(1))),
	}}
	var reference [][]*big.Int
	primes := fastStoragePrimes()
	for width := 1; width <= 3; width++ {
		input := newFastStorageCiphertextFromBig(t, params, values, level, width, false, nil)
		output, err := FastStorageRescale(input)
		require.NoError(t, err)
		want := roundBigIntMatrix(values, params.N(), q)
		requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
		require.Equal(t, width, output.StorageWidth())
		if reference == nil {
			reference = decodeFastStoragePolynomials(t, output)
		} else {
			requireBigIntMatricesEqual(t, reference, decodeFastStoragePolynomials(t, output))
		}
		for _, f := range primes[:width] {
			require.NotEqual(t, q, f, "q_ell must be independent of active F rows")
		}
	}
}

func TestFastStorageRescaleCoefficientAndNTTExactness(t *testing.T) {
	params := storageBoundaryParameters(t)
	level := 3
	q := params.Q()[level]
	qBig := new(big.Int).SetUint64(q)
	values := [][]*big.Int{
		{new(big.Int).Add(new(big.Int).Mul(qBig, big.NewInt(5)), big.NewInt(7)), new(big.Int).Neg(new(big.Int).Add(new(big.Int).Mul(qBig, big.NewInt(3)), big.NewInt(2))), big.NewInt(9)},
		{new(big.Int).Neg(new(big.Int).Add(new(big.Int).Mul(qBig, big.NewInt(2)), big.NewInt(4))), new(big.Int).Add(new(big.Int).Mul(qBig, big.NewInt(4)), big.NewInt(1))},
	}
	want := roundBigIntMatrix(values, params.N(), q)
	for _, ntt := range []bool{false, true} {
		input := newFastStorageCiphertextFromBig(t, params, values, level, 2, ntt, nil)
		output, err := FastStorageRescale(input)
		require.NoError(t, err)
		requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
		require.Equal(t, ntt, output.IsNTT())
		require.False(t, output.IsMontgomery())
	}
}

func TestFastStorageRescaleBoundsArePerComponentAndExact(t *testing.T) {
	params := storageBoundaryParameters(t)
	level := 2
	q := params.Q()[level]
	bounds := []*big.Int{
		new(big.Int).Add(new(big.Int).Mul(new(big.Int).SetUint64(q), big.NewInt(5)), big.NewInt(3)),
		new(big.Int).Add(new(big.Int).Mul(new(big.Int).SetUint64(q), big.NewInt(17)), big.NewInt(7)),
	}
	input := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{big.NewInt(1)}, {big.NewInt(-2)}}, level, 1, false, bounds)
	output, err := FastStorageRescale(input)
	require.NoError(t, err)
	want := make([]*big.Int, len(bounds))
	for component, bound := range bounds {
		want[component] = rescaleComponentBound(bound, q)
	}
	require.Equal(t, want, output.ComponentBounds())
	require.NotEqual(t, output.ComponentBounds()[0], output.ComponentBounds()[1], "nonuniform input bounds must remain per-component")
}

func TestFastStorageRescaleMetadataAndNoImplicitContraction(t *testing.T) {
	params := storageBoundaryParameters(t)
	level := 3
	q := params.Q()[level]
	values := [][]*big.Int{{big.NewInt(5)}, {big.NewInt(-9)}}
	for width := 1; width <= 3; width++ {
		input := newFastStorageCiphertextFromBig(t, params, values, level, width, width&1 == 1, nil)
		input.metadata.Scale = rlwe.NewScale(1 << 60)
		input.metadata.IsBatched = false
		input.metadata.IsBitReversed = true
		input.metadata.LogDimensions = ring.Dimensions{Rows: 1, Cols: 2}
		beforeMeta := input.MetaData()
		beforeParams := input.Parameters()
		output, err := FastStorageRescale(input)
		require.NoError(t, err)
		require.Equal(t, level-1, output.LogicalLevel())
		require.True(t, output.Scale().Equal(beforeMeta.Scale.Div(rlwe.NewScale(q))))
		require.Equal(t, input.Degree(), output.Degree())
		require.Equal(t, width, output.StorageWidth(), "logical level change must not contract private storage")
		require.Equal(t, input.IsNTT(), output.IsNTT())
		require.False(t, output.IsMontgomery())
		metadata := output.MetaData()
		require.Equal(t, beforeMeta.IsBatched, metadata.IsBatched)
		require.Equal(t, beforeMeta.IsBitReversed, metadata.IsBitReversed)
		require.Equal(t, beforeMeta.LogDimensions, metadata.LogDimensions)
		outputParams := output.Parameters()
		require.True(t, beforeParams.Equal(&outputParams))
	}

	bounds := []*big.Int{new(big.Int).Lsh(big.NewInt(1), 100)}
	wideInput := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{}}, level, 3, false, bounds)
	wideOutput, err := FastStorageRescale(wideInput)
	require.NoError(t, err)
	require.Equal(t, 3, wideOutput.StorageWidth())
	wantBound := rescaleComponentBound(bounds[0], q)
	capacity1, err := wideOutput.basis.centeredCapacity(1)
	require.NoError(t, err)
	require.Less(t, wantBound.Cmp(capacity1), 0, "post-Rescale bound would fit width 1, but width must remain unchanged")
}

func TestFastStorageRescaleStandardLogicalOracle(t *testing.T) {
	params := storageBoundaryParameters(t)
	for _, level := range []int{1, 3} {
		q := params.Q()[level]
		values := [][]*big.Int{
			{new(big.Int).Add(new(big.Int).Mul(new(big.Int).SetUint64(q), big.NewInt(2)), big.NewInt(3)), new(big.Int).Neg(new(big.Int).Add(new(big.Int).SetUint64(q), big.NewInt(1)))},
			{new(big.Int).Neg(new(big.Int).Add(new(big.Int).Mul(new(big.Int).SetUint64(q), big.NewInt(4)), big.NewInt(5)))},
		}
		for _, ntt := range []bool{false, true} {
			input := newFastStorageCiphertextFromBig(t, params, values, level, 3, ntt, nil)
			output, err := FastStorageRescale(input)
			require.NoError(t, err)
			assertStandardLogicalRescaleOracle(t, input, output)
			requireBigIntMatricesEqual(t, roundBigIntMatrix(values, params.N(), q), decodeFastStoragePolynomials(t, output))
		}
	}
}

func TestFastStorageRescaleNoncanonicalLiftCongruenceTheorem(t *testing.T) {
	params := storageBoundaryParameters(t)
	level := 1
	qProduct := big.NewInt(1)
	for i := 0; i <= level; i++ {
		qProduct.Mul(qProduct, new(big.Int).SetUint64(params.Q()[i]))
	}
	cValues := [][]int64{{3, -4, 7}, {-2, 5}}
	values := make([][]*big.Int, len(cValues))
	for component := range cValues {
		values[component] = make([]*big.Int, len(cValues[component]))
		for coefficient, c := range cValues[component] {
			k := int64(2)
			if (component+coefficient)&1 == 1 {
				k = -3
			}
			values[component][coefficient] = new(big.Int).Add(big.NewInt(c), new(big.Int).Mul(qProduct, big.NewInt(k)))
		}
	}
	q := params.Q()[level]
	for _, ntt := range []bool{false, true} {
		input := newFastStorageCiphertextFromBig(t, params, values, level, 3, ntt, nil)
		output, err := FastStorageRescale(input)
		require.NoError(t, err)
		assertStandardLogicalRescaleOracle(t, input, output)
		requireBigIntMatricesEqual(t, roundBigIntMatrix(values, params.N(), q), decodeFastStoragePolynomials(t, output))
	}
}

func TestFastStorageRescaleFailuresAreTransactional(t *testing.T) {
	params := storageBoundaryParameters(t)
	makeInput := func() *FastCiphertext {
		return newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{big.NewInt(7)}}, 2, 2, true, nil)
	}
	tests := []struct {
		name   string
		mutate func(*FastCiphertext)
		want   string
	}{
		{name: "level-zero", mutate: func(ct *FastCiphertext) { ct.logicalLevel = 0 }, want: "Level >= 1"},
		{name: "montgomery", mutate: func(ct *FastCiphertext) { ct.metadata.IsMontgomery = true }, want: "Montgomery"},
		{name: "malformed-bound", mutate: func(ct *FastCiphertext) { ct.componentBounds[0] = nil }, want: "bound cannot be nil"},
		{name: "malformed-row", mutate: func(ct *FastCiphertext) { ct.value[0].Coeffs[0] = ct.value[0].Coeffs[0][:1] }, want: "expected N"},
		{name: "non-positive-scale", mutate: func(ct *FastCiphertext) { ct.metadata.Scale = rlwe.NewScale(0) }, want: "scale must be positive"},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			input := makeInput()
			test.mutate(input)
			before := snapshotFastStorageCiphertext(input)
			output, err := FastStorageRescale(input)
			require.Nil(t, output)
			require.ErrorContains(t, err, test.want)
			require.Equal(t, before, snapshotFastStorageCiphertext(input))
		})
	}
	output, err := FastStorageRescale(nil)
	require.Nil(t, output)
	require.ErrorContains(t, err, "cannot be nil")
}

func TestFastStorageMulThenRescaleBoundAndOracle(t *testing.T) {
	params := storageBoundaryParameters(t)
	leftValues := [][]*big.Int{{big.NewInt(1), big.NewInt(-2)}, {big.NewInt(3)}}
	rightValues := [][]*big.Int{{big.NewInt(-1), big.NewInt(4)}, {big.NewInt(2), big.NewInt(-3)}}
	left := newFastStorageCiphertextFromBig(t, params, leftValues, 2, 1, true, nil)
	right := newFastStorageCiphertextFromBig(t, params, rightValues, 2, 1, true, nil)
	left.metadata.Scale = rlwe.NewScale(1 << 60)
	right.metadata.Scale = rlwe.NewScale(1 << 60)
	mul, err := FastStorageMul(left, right)
	require.NoError(t, err)
	inputBounds := mul.ComponentBounds()
	rescaled, err := FastStorageRescale(mul)
	require.NoError(t, err)
	q := params.Q()[2]
	wantBounds := make([]*big.Int, len(inputBounds))
	for component, bound := range inputBounds {
		wantBounds[component] = rescaleComponentBound(bound, q)
	}
	gotBounds := rescaled.ComponentBounds()
	require.Len(t, gotBounds, len(wantBounds))
	for component := range wantBounds {
		require.Zero(t, wantBounds[component].Cmp(gotBounds[component]), "Rescale must consume Mul provenance without rescanning, component=%d", component)
	}
	require.Equal(t, 2, rescaled.Degree(), "raw multiplication degree must remain unrelinearized")
	require.Equal(t, 1, rescaled.LogicalLevel())
	require.True(t, rescaled.Scale().Equal(mul.Scale().Div(rlwe.NewScale(q))))

	product := expectedNegacyclicBigCiphertextProduct(leftValues, rightValues, params.N())
	requireBigIntMatricesEqual(t, roundBigIntMatrix(product, params.N(), q), decodeFastStoragePolynomials(t, rescaled))
	assertStandardLogicalRescaleOracle(t, mul, rescaled)
}

func newFastStorageCiphertextFromBig(t *testing.T, params ckks.Parameters, polynomials [][]*big.Int, level, width int, ntt bool, bounds []*big.Int) *FastCiphertext {
	t.Helper()
	require.NotEmpty(t, polynomials)
	ct, err := NewFastCiphertext(params, len(polynomials)-1, level, width)
	require.NoError(t, err)
	ct.metadata.IsNTT = ntt
	ct.componentBounds = zeroComponentBounds(len(polynomials))
	for component, polynomial := range polynomials {
		for coefficient, value := range polynomial {
			require.NotNil(t, value)
			require.Less(t, coefficient, params.N())
			integer := storageIntegerFromBig(value)
			encoded, err := encodeCenteredStorageValue(ct.basis, integer, width)
			require.NoError(t, err)
			for row := 0; row < width; row++ {
				ct.value[component].Coeffs[row][coefficient] = encoded[row]
			}
			magnitude := new(big.Int).Abs(new(big.Int).Set(value))
			if magnitude.Cmp(ct.componentBounds[component]) > 0 {
				ct.componentBounds[component].Set(magnitude)
			}
		}
	}
	if bounds != nil {
		ct.componentBounds = cloneComponentBounds(bounds)
	}
	if ntt {
		for component := range ct.value {
			for row := 0; row < width; row++ {
				subring, err := ct.basis.subring(row)
				require.NoError(t, err)
				subring.NTT(ct.value[component].Coeffs[row], ct.value[component].Coeffs[row])
			}
		}
	}
	require.NoError(t, validateFastStorageState(ct))
	return ct
}

func storageIntegerFromBig(value *big.Int) storageInteger {
	magnitude := new(big.Int).Abs(new(big.Int).Set(value))
	return storageInteger{magnitude: storageFromBig(magnitude), negative: value.Sign() < 0}
}

func roundBigIntByUint64(value *big.Int, divisor uint64) *big.Int {
	if value.Sign() == 0 {
		return new(big.Int)
	}
	magnitude := new(big.Int).Abs(new(big.Int).Set(value))
	q := new(big.Int).SetUint64(divisor)
	magnitude.Add(magnitude, new(big.Int).SetUint64((divisor-1)/2))
	magnitude.Quo(magnitude, q)
	if value.Sign() < 0 {
		magnitude.Neg(magnitude)
	}
	return magnitude
}

func roundBigIntMatrix(polynomials [][]*big.Int, n int, divisor uint64) [][]*big.Int {
	want := make([][]*big.Int, len(polynomials))
	for component, polynomial := range polynomials {
		want[component] = make([]*big.Int, n)
		for coefficient := 0; coefficient < n; coefficient++ {
			value := new(big.Int)
			if coefficient < len(polynomial) {
				value.Set(polynomial[coefficient])
			}
			want[component][coefficient] = roundBigIntByUint64(value, divisor)
		}
	}
	return want
}

func expectedNegacyclicBigCiphertextProduct(left, right [][]*big.Int, n int) [][]*big.Int {
	want := make([][]*big.Int, len(left)+len(right)-1)
	for component := range want {
		want[component] = make([]*big.Int, n)
		for coefficient := range want[component] {
			want[component][coefficient] = new(big.Int)
		}
	}
	for leftComponent := range left {
		for rightComponent := range right {
			for i, a := range left[leftComponent] {
				for j, b := range right[rightComponent] {
					index := i + j
					term := new(big.Int).Mul(a, b)
					if index >= n {
						index -= n
						term.Neg(term)
					}
					want[leftComponent+rightComponent][index].Add(want[leftComponent+rightComponent][index], term)
				}
			}
		}
	}
	return want
}

func assertStandardLogicalRescaleOracle(t *testing.T, input, fastOutput *FastCiphertext) {
	t.Helper()
	logicalInput, err := input.ExportToLogical(FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	standardOutput := ckks.NewCiphertext(input.params, input.Degree(), input.LogicalLevel()-1)
	standardEvaluator := ckks.NewEvaluator(input.params, nil)
	require.NoError(t, standardEvaluator.Rescale(logicalInput, standardOutput))
	logicalFastOutput, err := fastOutput.ExportToLogical(FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	require.Equal(t, input.LogicalLevel()-1, standardOutput.Level())
	require.Equal(t, standardOutput.Level(), logicalFastOutput.Level())
	require.True(t, standardOutput.Scale.Equal(logicalFastOutput.Scale))
	require.Equal(t, standardOutput.Degree(), logicalFastOutput.Degree())
	for component := range standardOutput.Value {
		for level := 0; level <= standardOutput.Level(); level++ {
			require.Equal(t, standardOutput.Value[component].Coeffs[level], logicalFastOutput.Value[component].Coeffs[level], "component=%d level=%d", component, level)
		}
	}
}
