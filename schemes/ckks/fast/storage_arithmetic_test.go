package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFastStorageRequiredWidthStrictExactThresholds(t *testing.T) {
	params := storageBoundaryParameters(t)
	basis, err := newFastStorageBasis(params.LogN())
	require.NoError(t, err)

	capacity1, err := basis.centeredCapacity(1)
	require.NoError(t, err)
	width, err := requiredFastStorageWidth(basis, []*big.Int{capacity1})
	require.NoError(t, err)
	require.Equal(t, 1, width)

	product1, err := basis.product(1)
	require.NoError(t, err)
	strictFailure := new(big.Int).Add(new(big.Int).Rsh(product1, 1), big.NewInt(1))
	width, err = requiredFastStorageWidth(basis, []*big.Int{strictFailure})
	require.NoError(t, err)
	require.Equal(t, 2, width, "2B >= S1 must not be accepted at width 1")

	capacity2, err := basis.centeredCapacity(2)
	require.NoError(t, err)
	width, err = requiredFastStorageWidth(basis, []*big.Int{new(big.Int).Add(capacity2, big.NewInt(1))})
	require.NoError(t, err)
	require.Equal(t, 3, width)

	capacity3, err := basis.centeredCapacity(3)
	require.NoError(t, err)
	_, err = requiredFastStorageWidth(basis, []*big.Int{new(big.Int).Add(capacity3, big.NewInt(1))})
	require.ErrorContains(t, err, "all Fast storage widths")
	_, err = requiredFastStorageWidth(basis, []*big.Int{big.NewInt(-1)})
	require.ErrorContains(t, err, "cannot be negative")
}

func TestFastStorageExpansionPreservesExactLiftAllWideningsAndDomains(t *testing.T) {
	params := storageBoundaryParameters(t)
	polynomials := [][]int64{{0, 1, -1, 17, -23}, {8, -5, 2}}
	for _, test := range []struct {
		from, to int
	}{{1, 2}, {1, 3}, {2, 3}} {
		for _, ntt := range []bool{false, true} {
			name := "coeff"
			if ntt {
				name = "ntt"
			}
			t.Run(name+"/"+string(rune('0'+test.from))+"-to-"+string(rune('0'+test.to)), func(t *testing.T) {
				input := newSyntheticStorageCiphertext(t, params, polynomials, 3, test.from, ntt, nil)
				input.metadata.Scale = rlwe.NewScale(123)
				want := decodeFastStoragePolynomials(t, input)
				output, err := input.ExpandStorage(test.to)
				require.NoError(t, err)
				require.Equal(t, test.to, output.StorageWidth())
				require.Equal(t, input.LogicalLevel(), output.LogicalLevel())
				require.Equal(t, input.Degree(), output.Degree())
				require.True(t, input.Scale().Equal(output.Scale()))
				require.Equal(t, input.IsNTT(), output.IsNTT())
				require.Equal(t, input.ComponentBounds(), output.ComponentBounds())
				require.Equal(t, want, decodeFastStoragePolynomials(t, output))
				require.NoError(t, validateFastStorageState(output))
				require.Equal(t, test.from, input.StorageWidth(), "expansion must not mutate the input")
			})
		}
	}
}

func TestFastStorageAddSubExactnessMetadataAndMixedWidths(t *testing.T) {
	params := storageBoundaryParameters(t)
	leftPolys := [][]int64{{3, -1, 2}}
	rightPolys := [][]int64{{-2, 4}, {5, -3, 1}}
	for _, ntt := range []bool{false, true} {
		inputDomain := "coeff"
		if ntt {
			inputDomain = "ntt"
		}
		t.Run(inputDomain, func(t *testing.T) {
			left := newSyntheticStorageCiphertext(t, params, leftPolys, 4, 1, ntt, nil)
			right := newSyntheticStorageCiphertext(t, params, rightPolys, 2, 2, ntt, nil)
			left.metadata.Scale = rlwe.NewScale(101)
			right.metadata.Scale = rlwe.NewScale(101)
			leftBefore, rightBefore := snapshotFastStorageCiphertext(left), snapshotFastStorageCiphertext(right)

			for _, subtract := range []bool{false, true} {
				var got *FastCiphertext
				var err error
				if subtract {
					got, err = FastStorageSub(left, right)
				} else {
					got, err = FastStorageAdd(left, right)
				}
				require.NoError(t, err)
				require.Equal(t, 1, got.Degree())
				require.Equal(t, 2, got.LogicalLevel())
				require.Equal(t, 2, got.StorageWidth(), "mixed inputs must expand to the wider input basis")
				require.Equal(t, ntt, got.IsNTT())
				require.True(t, got.Scale().Equal(left.Scale()))
				want := expectedAddSubPolynomials(leftPolys, rightPolys, params.N(), subtract)
				requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, got))
				require.Equal(t, expectedSummedBounds(left.ComponentBounds(), right.ComponentBounds(), 2), got.ComponentBounds())
				assertLogicalExportMatches(t, got, want)
			}
			require.Equal(t, leftBefore, snapshotFastStorageCiphertext(left))
			require.Equal(t, rightBefore, snapshotFastStorageCiphertext(right))
		})
	}

	left := newSyntheticStorageCiphertext(t, params, leftPolys, 0, 1, false, nil)
	right := newSyntheticStorageCiphertext(t, params, leftPolys, 0, 1, false, nil)
	right.metadata.Scale = rlwe.NewScale(102)
	_, err := FastStorageAdd(left, right)
	require.ErrorContains(t, err, "equal CKKS scales")
	right.metadata.Scale = left.metadata.Scale
	right.metadata.IsBatched = !left.metadata.IsBatched
	_, err = FastStorageAdd(left, right)
	require.ErrorContains(t, err, "compatible plaintext metadata")
}

func TestFastStorageRawMulExactnessBoundsScaleAndLogicalExport(t *testing.T) {
	params := storageBoundaryParameters(t)
	leftPolys := [][]int64{{1, 2, -1}, {2, -1}}
	rightPolys := [][]int64{{3, -2}, {1, 4, 2}}
	left := newSyntheticStorageCiphertext(t, params, leftPolys, 3, 1, true, nil)
	right := newSyntheticStorageCiphertext(t, params, rightPolys, 2, 1, true, nil)
	left.metadata.Scale = rlwe.NewScale(5)
	right.metadata.Scale = rlwe.NewScale(7)
	leftBefore, rightBefore := snapshotFastStorageCiphertext(left), snapshotFastStorageCiphertext(right)

	got, err := FastStorageMul(left, right)
	require.NoError(t, err)
	require.Equal(t, 2, got.Degree())
	require.Equal(t, 2, got.LogicalLevel())
	require.Equal(t, 1, got.StorageWidth())
	require.True(t, got.IsNTT())
	require.True(t, got.Scale().Equal(rlwe.NewScale(35)))
	want := expectedNegacyclicCiphertextProduct(leftPolys, rightPolys, params.N())
	requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, got))

	leftBounds, rightBounds := left.ComponentBounds(), right.ComponentBounds()
	wantBounds := []*big.Int{
		new(big.Int).Mul(leftBounds[0], rightBounds[0]),
		new(big.Int).Add(new(big.Int).Mul(leftBounds[0], rightBounds[1]), new(big.Int).Mul(leftBounds[1], rightBounds[0])),
		new(big.Int).Mul(leftBounds[1], rightBounds[1]),
	}
	for _, bound := range wantBounds {
		bound.Mul(bound, new(big.Int).SetUint64(uint64(params.N())))
	}
	require.Equal(t, wantBounds, got.ComponentBounds(), "per-component nonuniform convolution bounds")
	assertLogicalExportMatches(t, got, want)
	require.Equal(t, leftBefore, snapshotFastStorageCiphertext(left))
	require.Equal(t, rightBefore, snapshotFastStorageCiphertext(right))
}

func TestFastStorageMulPlansWidthsWithoutContraction(t *testing.T) {
	params := storageBoundaryParameters(t)
	cases := []struct {
		name       string
		inputWidth int
		bound      *big.Int
		wantWidth  int
	}{
		{name: "width-1", inputWidth: 1, bound: big.NewInt(1), wantWidth: 1},
		{name: "expand-to-width-2", inputWidth: 1, bound: new(big.Int).Lsh(big.NewInt(1), 30), wantWidth: 2},
		{name: "expand-to-width-3", inputWidth: 2, bound: new(big.Int).Lsh(big.NewInt(1), 60), wantWidth: 3},
		{name: "no-contraction", inputWidth: 3, bound: big.NewInt(1), wantWidth: 3},
	}
	for _, test := range cases {
		t.Run(test.name, func(t *testing.T) {
			bounds := []*big.Int{new(big.Int).Set(test.bound)}
			left := newSyntheticStorageCiphertext(t, params, [][]int64{{}}, 0, test.inputWidth, true, bounds)
			right := newSyntheticStorageCiphertext(t, params, [][]int64{{}}, 0, test.inputWidth, true, bounds)
			got, err := FastStorageMul(left, right)
			require.NoError(t, err)
			require.Equal(t, test.wantWidth, got.StorageWidth())
			require.NoError(t, validateFastStorageState(got))
		})
	}
}

func TestFastStorageArithmeticCapacityFailureIsTransactional(t *testing.T) {
	params := storageBoundaryParameters(t)
	basis, err := newFastStorageBasis(params.LogN())
	require.NoError(t, err)
	capacity, err := basis.centeredCapacity(3)
	require.NoError(t, err)
	bound := new(big.Int).Mul(capacity, big.NewInt(2))
	bound.Quo(bound, big.NewInt(3))
	override := []*big.Int{bound}
	left := newSyntheticStorageCiphertext(t, params, [][]int64{{}}, 1, 3, false, override)
	right := newSyntheticStorageCiphertext(t, params, [][]int64{{}}, 1, 3, false, override)
	leftBefore, rightBefore := snapshotFastStorageCiphertext(left), snapshotFastStorageCiphertext(right)

	got, err := FastStorageAdd(left, right)
	require.Nil(t, got)
	require.ErrorContains(t, err, "all Fast storage widths")
	require.Equal(t, leftBefore, snapshotFastStorageCiphertext(left))
	require.Equal(t, rightBefore, snapshotFastStorageCiphertext(right))
}

func TestFastStorageArithmeticRejectsInvalidDomainsParametersAndMetadata(t *testing.T) {
	params := storageBoundaryParameters(t)
	coeff := newSyntheticStorageCiphertext(t, params, [][]int64{{1}}, 0, 1, false, nil)
	ntt := newSyntheticStorageCiphertext(t, params, [][]int64{{1}}, 0, 1, true, nil)
	_, err := FastStorageMul(coeff, coeff)
	require.ErrorContains(t, err, "requires NTT-domain")
	_, err = FastStorageAdd(coeff, ntt)
	require.ErrorContains(t, err, "matching coefficient/NTT domains")

	montgomery := coeff.CopyNew()
	montgomery.metadata.IsMontgomery = true
	_, err = FastStorageAdd(montgomery, coeff)
	require.ErrorContains(t, err, "Montgomery")

	otherParams := storageBoundaryParametersLogN(t, 5)
	other := newSyntheticStorageCiphertext(t, otherParams, [][]int64{{1}}, 0, 1, false, nil)
	_, err = FastStorageAdd(coeff, other)
	require.ErrorContains(t, err, "matching CKKS parameters")

	incompatible := coeff.CopyNew()
	incompatible.metadata.LogDimensions.Rows++
	_, err = FastStorageAdd(coeff, incompatible)
	require.ErrorContains(t, err, "compatible plaintext metadata")

	malformedRows := coeff.CopyNew()
	malformedRows.value[0].Coeffs[0] = malformedRows.value[0].Coeffs[0][:1]
	_, err = FastStorageAdd(malformedRows, coeff)
	require.ErrorContains(t, err, "expected N")

	malformedBounds := coeff.CopyNew()
	malformedBounds.componentBounds[0] = nil
	_, err = FastStorageAdd(malformedBounds, coeff)
	require.ErrorContains(t, err, "bound cannot be nil")

	malformedBounds = coeff.CopyNew()
	malformedBounds.componentBounds[0].Neg(big.NewInt(1))
	_, err = FastStorageAdd(malformedBounds, coeff)
	require.ErrorContains(t, err, "bound cannot be negative")
}

func TestFastStorageImportRecordsExactCoefficientAndNTTBounds(t *testing.T) {
	params := storageBoundaryParameters(t)
	source, expected := storageBoundarySource(t, params)
	for _, useNTT := range []bool{false, true} {
		input := source.CopyNew()
		if useNTT {
			logicalRing := params.RingQ().AtLevel(0)
			for component := range input.Value {
				logicalRing.SubRings[0].NTT(input.Value[component].Coeffs[0], input.Value[component].Coeffs[0])
			}
			input.IsNTT = true
		}
		fastCT, err := ImportLevel0(params, input, 2, FastCiphertextDomain{IsNTT: useNTT})
		require.NoError(t, err)
		for component := range expected {
			want := maxAbsoluteBigInt(expected[component])
			require.Equal(t, want, fastCT.ComponentBounds()[component])
		}
	}
}

func newSyntheticStorageCiphertext(t *testing.T, params ckks.Parameters, polynomials [][]int64, level, width int, ntt bool, bounds []*big.Int) *FastCiphertext {
	t.Helper()
	ct, err := NewFastCiphertext(params, len(polynomials)-1, level, width)
	require.NoError(t, err)
	ct.metadata.IsNTT = ntt
	ct.metadata.IsMontgomery = false
	ct.componentBounds = zeroComponentBounds(len(polynomials))
	for component, polynomial := range polynomials {
		for coefficient, value := range polynomial {
			require.Less(t, coefficient, params.N())
			integer := big.NewInt(value)
			encoded, err := ct.basis.encode(integer, width)
			require.NoError(t, err)
			for row := 0; row < width; row++ {
				ct.value[component].Coeffs[row][coefficient] = encoded[row]
			}
			if new(big.Int).Abs(new(big.Int).Set(integer)).Cmp(ct.componentBounds[component]) > 0 {
				ct.componentBounds[component].Abs(integer)
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

func decodeFastStoragePolynomials(t *testing.T, ct *FastCiphertext) [][]*big.Int {
	t.Helper()
	decoded := make([][]*big.Int, len(ct.value))
	for component := range ct.value {
		coeffRows := make([][]uint64, ct.StorageWidth())
		for row := 0; row < ct.StorageWidth(); row++ {
			coeffRows[row] = make([]uint64, ct.N())
			subring, err := ct.basis.subring(row)
			require.NoError(t, err)
			if ct.IsNTT() {
				subring.INTT(ct.value[component].Coeffs[row], coeffRows[row])
			} else {
				copy(coeffRows[row], ct.value[component].Coeffs[row])
			}
		}
		decoded[component] = make([]*big.Int, ct.N())
		for coefficient := 0; coefficient < ct.N(); coefficient++ {
			var residues [3]uint64
			for row := 0; row < ct.StorageWidth(); row++ {
				residues[row] = coeffRows[row][coefficient]
			}
			value, err := ct.basis.decodeFixed(residues, ct.StorageWidth())
			require.NoError(t, err, "component=%d coefficient=%d residues=%v", component, coefficient, residues)
			decoded[component][coefficient] = signedStorageBig(value)
		}
	}
	return decoded
}

func expectedAddSubPolynomials(left, right [][]int64, n int, subtract bool) [][]*big.Int {
	degree := max(len(left), len(right))
	want := make([][]*big.Int, degree)
	for component := 0; component < degree; component++ {
		want[component] = make([]*big.Int, n)
		for coefficient := 0; coefficient < n; coefficient++ {
			value := big.NewInt(0)
			if component < len(left) && coefficient < len(left[component]) {
				value.SetInt64(left[component][coefficient])
			}
			if component < len(right) && coefficient < len(right[component]) {
				if subtract {
					value.Sub(value, big.NewInt(right[component][coefficient]))
				} else {
					value.Add(value, big.NewInt(right[component][coefficient]))
				}
			}
			want[component][coefficient] = value
		}
	}
	return want
}

func expectedNegacyclicCiphertextProduct(left, right [][]int64, n int) [][]*big.Int {
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
					term := big.NewInt(a * b)
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

func expectedSummedBounds(left, right []*big.Int, count int) []*big.Int {
	want := zeroComponentBounds(count)
	for component := range want {
		if component < len(left) {
			want[component].Add(want[component], left[component])
		}
		if component < len(right) {
			want[component].Add(want[component], right[component])
		}
	}
	return want
}

func maxAbsoluteBigInt(values []*big.Int) *big.Int {
	maximum := new(big.Int)
	for _, value := range values {
		magnitude := new(big.Int).Abs(new(big.Int).Set(value))
		if magnitude.Cmp(maximum) > 0 {
			maximum.Set(magnitude)
		}
	}
	return maximum
}

func requireBigIntMatricesEqual(t *testing.T, expected, actual [][]*big.Int) {
	t.Helper()
	require.Len(t, actual, len(expected))
	for component := range expected {
		require.Len(t, actual[component], len(expected[component]))
		for coefficient := range expected[component] {
			require.Zero(t, expected[component][coefficient].Cmp(actual[component][coefficient]), "component=%d coefficient=%d expected=%s actual=%s", component, coefficient, expected[component][coefficient], actual[component][coefficient])
		}
	}
}

type fastStorageSnapshot struct {
	rows   []ring.Poly
	level  int
	width  int
	bounds []*big.Int
	meta   rlwe.MetaData
}

func snapshotFastStorageCiphertext(ct *FastCiphertext) fastStorageSnapshot {
	return fastStorageSnapshot{
		rows:   cloneStorageRows(ct.value),
		level:  ct.logicalLevel,
		width:  ct.activeStorageWidth,
		bounds: cloneComponentBounds(ct.componentBounds),
		meta:   cloneFastMetadata(ct.metadata),
	}
}

func assertLogicalExportMatches(t *testing.T, ct *FastCiphertext, expected [][]*big.Int) {
	t.Helper()
	logical, err := ct.ExportToLogical(FastCiphertextDomain{})
	require.NoError(t, err)
	for component := range expected {
		for coefficient, value := range expected[component] {
			for level := 0; level <= ct.LogicalLevel(); level++ {
				modulus := new(big.Int).SetUint64(ct.params.Q()[level])
				want := new(big.Int).Mod(new(big.Int).Set(value), modulus).Uint64()
				require.Equal(t, want, logical.Value[component].Coeffs[level][coefficient])
			}
		}
	}
}

func storageBoundaryParametersLogN(t *testing.T, logN int) ckks.Parameters {
	t.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(45, 2<<uint(logN))
	q, err := generator.NextAlternatingPrimes(3)
	require.NoError(t, err)
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: logN, Q: q, LogDefaultScale: 30})
	require.NoError(t, err)
	return params
}
