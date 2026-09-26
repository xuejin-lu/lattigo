package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFastStoragePlaintextMirrorRoundTripAndBounds(t *testing.T) {
	params := fastStoragePlaintextTestParameters(t)
	level := params.MaxLevel()
	basis, err := fastStorageBasisForLogN(params.LogN())
	require.NoError(t, err)
	capacity, err := basis.centeredCapacity(3)
	require.NoError(t, err)
	q0 := new(big.Int).SetUint64(params.Q()[0])
	values := make([]*big.Int, params.N())
	values[0] = big.NewInt(0)
	values[1] = big.NewInt(1)
	values[2] = big.NewInt(-1)
	values[3] = new(big.Int).Add(q0, big.NewInt(19))
	values[4] = new(big.Int).Neg(new(big.Int).Add(q0, big.NewInt(23)))
	values[5] = new(big.Int).Sub(capacity, big.NewInt(17))
	values[6] = new(big.Int).Neg(new(big.Int).Set(values[5]))
	for i := 7; i < len(values); i++ {
		values[i] = new(big.Int)
	}

	logicalQ := logicalQPolynomialFromIntegers(params, level, values)
	metadata := rlwe.PlaintextMetaData{Scale: rlwe.NewScale(1 << 20), IsBatched: true, LogDimensions: ring.Dimensions{Cols: 3}}
	pt, err := NewFastStoragePlaintextMirror(params, logicalQ, level, metadata, true, true)
	require.NoError(t, err)
	require.Equal(t, 3, pt.StorageWidth())
	require.True(t, pt.IsNTT())
	require.False(t, pt.IsMontgomery())
	require.Equal(t, level, pt.LevelQ())
	require.True(t, metadata.Scale.Equal(pt.Scale()))
	require.Equal(t, values[5], pt.MaxAbsCoefficient())
	wantL1 := new(big.Int)
	for _, value := range values {
		wantL1.Add(wantL1, new(big.Int).Abs(new(big.Int).Set(value)))
	}
	require.Equal(t, wantL1, pt.L1Norm())
	for row := 0; row < 3; row++ {
		coefficients := make([]uint64, params.N())
		subring, _ := basis.subring(row)
		subring.INTT(pt.value.Coeffs[row], coefficients)
		for coefficient, value := range values {
			want := new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(fastStoragePrimes()[row])).Uint64()
			require.Equal(t, want, coefficients[coefficient], "F row %d coefficient %d", row, coefficient)
		}
	}

	tooLarge := append([]*big.Int(nil), values...)
	tooLarge[5] = new(big.Int).Add(capacity, big.NewInt(1))
	_, err = NewFastStoragePlaintextMirror(params, logicalQPolynomialFromIntegers(params, level, tooLarge), level, metadata, true, true)
	require.ErrorContains(t, err, "strict width-3 centered capacity")

	incomplete := ring.NewPoly(params.N(), level-1)
	_, err = NewFastStoragePlaintextMirror(params, incomplete, level, metadata, true, true)
	require.ErrorContains(t, err, "expected")

	for _, profile := range []struct {
		name  string
		level int
	}{
		{name: "q01", level: 1},
		{name: "q012", level: 2},
	} {
		t.Run(profile.name, func(t *testing.T) {
			profileValues := make([]*big.Int, params.N())
			for i := range profileValues {
				profileValues[i] = new(big.Int)
			}
			if profile.name == "q01" {
				profileValues[0].Add(q0, big.NewInt(31))
				profileValues[1].Neg(new(big.Int).Add(q0, big.NewInt(9)))
			} else {
				q01 := new(big.Int).Mul(new(big.Int).SetUint64(params.Q()[0]), new(big.Int).SetUint64(params.Q()[1]))
				profileValues[0].Add(q01, big.NewInt(37))
				profileValues[1].Neg(new(big.Int).Add(q01, big.NewInt(13)))
			}
			profileMetadata := metadata
			profile, err := NewFastStoragePlaintextMirror(params, logicalQPolynomialFromIntegers(params, profile.level, profileValues), profile.level, profileMetadata, true, true)
			require.NoError(t, err)
			require.Equal(t, profileValues[0], decodeFastStoragePlaintext(t, profile)[0])
		})
	}
}

func TestFastStorageMulPlaintextMatchesIndependentNegacyclicConvolution(t *testing.T) {
	params := fastStoragePlaintextTestParameters(t)
	level := 2
	inputPolys := [][]*big.Int{{big.NewInt(3), big.NewInt(-2), big.NewInt(5)}, {big.NewInt(-7), big.NewInt(4)}}
	ct := newFastStorageCiphertextFromBig(t, params, inputPolys, level, 3, true, nil)
	ct.metadata.LogDimensions = ring.Dimensions{Cols: 3}
	plain := make([]*big.Int, params.N())
	for i := range plain {
		plain[i] = new(big.Int)
	}
	plain[0].SetInt64(2)
	plain[1].SetInt64(-3)
	plain[4].SetInt64(5)
	metadata := ct.metadata.PlaintextMetaData
	metadata.Scale = rlwe.NewScale(32)
	pt, err := NewFastStoragePlaintextMirror(params, logicalQPolynomialFromIntegers(params, level, plain), level, metadata, true, true)
	require.NoError(t, err)

	got, err := FastStorageMulPlaintext(ct, pt)
	require.NoError(t, err)
	want := make([][]*big.Int, ct.Degree()+1)
	for component := range inputPolys {
		want[component] = independentNegacyclicProduct(inputPolys[component], plain, params.N())
	}
	requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, got))
	for component, bound := range ct.ComponentBounds() {
		require.Equal(t, new(big.Int).Mul(bound, pt.L1Norm()), got.ComponentBounds()[component])
	}
	require.True(t, got.Scale().Equal(ct.Scale().Mul(pt.Scale())))
	require.Equal(t, ct.LogicalLevel(), got.LogicalLevel())
	require.Equal(t, 3, got.StorageWidth())
}

func TestFastStorageAutomorphismMatchesIndependentSignedPermutation(t *testing.T) {
	params := fastStoragePlaintextTestParameters(t)
	input := [][]*big.Int{{big.NewInt(2), big.NewInt(-3), big.NewInt(5)}, {big.NewInt(-7), big.NewInt(11)}}
	ct := newFastStorageCiphertextFromBig(t, params, input, 2, 3, true, nil)
	ct.metadata.Scale = rlwe.NewScale(1234)
	galEl := params.GaloisElement(3)
	got, err := FastStorageAutomorphism(ct, galEl)
	require.NoError(t, err)
	want := make([][]*big.Int, len(input))
	for component := range input {
		want[component] = independentAutomorphism(input[component], params.N(), galEl)
	}
	requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, got))
	require.Equal(t, ct.ComponentBounds(), got.ComponentBounds())
	require.Equal(t, ct.LogicalLevel(), got.LogicalLevel())
	require.True(t, ct.Scale().Equal(got.Scale()))
	require.Equal(t, ct.MetaData(), got.MetaData())
}

func TestFastStorageLinearTransformDirectAndBSGSExactness(t *testing.T) {
	params := fastStoragePlaintextTestParameters(t)
	level := 2
	input := [][]*big.Int{{big.NewInt(3), big.NewInt(-2), big.NewInt(5)}, {big.NewInt(-7), big.NewInt(4)}}
	ct := newFastStorageCiphertextFromBig(t, params, input, level, 3, true, nil)
	ct.metadata.LogDimensions = ring.Dimensions{Cols: 2}
	diagonals := map[int][]complex128{
		0: {2, -1, 0, 3},
		1: {1, 0, 2, -2},
		3: {-1, 2, 1, 0},
	}
	direct := encodedTestLinearTransformation(t, params, level, diagonals, -1)
	bsgs := encodedTestLinearTransformation(t, params, level, diagonals, 1)
	directMirror, err := NewFastStorageLinearTransformation(params, direct)
	require.NoError(t, err)
	bsgsMirror, err := NewFastStorageLinearTransformation(params, bsgs)
	require.NoError(t, err)
	require.Equal(t, len(direct.Vec), directMirror.DiagonalCount())
	require.Equal(t, len(bsgs.Vec), bsgsMirror.DiagonalCount())
	require.NotZero(t, bsgsMirror.N1())

	gotDirect, err := FastStorageLinearTransform(ct, directMirror)
	require.NoError(t, err)
	want := make([][]*big.Int, len(input))
	for component := range input {
		want[component] = make([]*big.Int, params.N())
		for i := range want[component] {
			want[component][i] = new(big.Int)
		}
	}
	for diagonal, plaintext := range directMirror.diagonals {
		plainCoefficients := decodeFastStoragePlaintext(t, plaintext)
		for component := range input {
			rotated := independentAutomorphism(input[component], params.N(), params.GaloisElement(diagonal))
			term := independentNegacyclicProduct(rotated, plainCoefficients, params.N())
			for coefficient := range term {
				want[component][coefficient].Add(want[component][coefficient], term[coefficient])
			}
		}
	}
	requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, gotDirect))

	gotBSGS, err := FastStorageLinearTransform(ct, bsgsMirror)
	require.NoError(t, err)
	requireBigIntMatricesEqual(t, decodeFastStoragePolynomials(t, gotDirect), decodeFastStoragePolynomials(t, gotBSGS))
	require.True(t, ct.Scale().Mul(direct.Scale).Equal(gotDirect.Scale()))
	require.Equal(t, ct.LogicalLevel(), gotDirect.LogicalLevel())
}

func TestFastStorageLinearTransformCapacityFailureIsTransactional(t *testing.T) {
	params := fastStoragePlaintextTestParameters(t)
	basis, err := fastStorageBasisForLogN(params.LogN())
	require.NoError(t, err)
	capacity, err := basis.centeredCapacity(3)
	require.NoError(t, err)
	inputValue := new(big.Int).Rsh(new(big.Int).Set(capacity), 1)
	ct := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{inputValue}, {new(big.Int)}}, 2, 3, true, []*big.Int{inputValue, new(big.Int)})
	ct.metadata.LogDimensions = ring.Dimensions{Cols: 1}
	before := snapshotFastStorageCiphertext(ct)
	plain := make([]*big.Int, params.N())
	for i := range plain {
		plain[i] = new(big.Int)
	}
	plain[0].SetInt64(5)
	lt := logicalTestLinearTransformation(t, params, 2, map[int][]*big.Int{0: plain}, ring.Dimensions{Cols: 1})
	matrix, err := NewFastStorageLinearTransformation(params, lt)
	require.NoError(t, err)
	require.ErrorContains(t, errIfCapacityFails(ct, matrix), "capacity preflight")
	require.Equal(t, inputValue, ct.ComponentBounds()[0])
	require.Equal(t, before, snapshotFastStorageCiphertext(ct))
	require.NoError(t, validateFastStorageState(ct))
}

func errIfCapacityFails(ct *FastCiphertext, matrix *FastStorageLinearTransformation) error {
	_, err := FastStorageLinearTransform(ct, matrix)
	return err
}

func fastStoragePlaintextTestParameters(t *testing.T) ckks.Parameters {
	t.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(50, 2*16)
	q, err := generator.NextAlternatingPrimes(5)
	require.NoError(t, err)
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 4, Q: q, LogP: []int{50}, LogDefaultScale: 30})
	require.NoError(t, err)
	return params
}

func logicalQPolynomialFromIntegers(params ckks.Parameters, level int, integers []*big.Int) ring.Poly {
	poly := ring.NewPoly(params.N(), level)
	for row := 0; row <= level; row++ {
		q := new(big.Int).SetUint64(params.Q()[row])
		for coefficient, integer := range integers {
			poly.Coeffs[row][coefficient] = new(big.Int).Mod(new(big.Int).Set(integer), q).Uint64()
		}
	}
	logicalRing := params.RingQ().AtLevel(level)
	logicalRing.NTT(poly, poly)
	logicalRing.MForm(poly, poly)
	return poly
}

func encodedTestLinearTransformation(t *testing.T, params ckks.Parameters, level int, diagonals map[int][]complex128, ratio int) lintrans.LinearTransformation {
	t.Helper()
	indices := make([]int, 0, len(diagonals))
	for diagonal := range diagonals {
		indices = append(indices, diagonal)
	}
	lt := lintrans.NewLinearTransformation(params, lintrans.Parameters{
		DiagonalsIndexList: indices, LevelQ: level, LevelP: params.MaxLevelP(),
		Scale: rlwe.NewScale(1 << 20), LogDimensions: ring.Dimensions{Cols: 2},
		LogBabyStepGiantStepRatio: ratio,
	})
	encoder := ckks.NewEncoder(params)
	if err := lintrans.Encode(encoder, lintrans.Diagonals[complex128](diagonals), lt); err != nil {
		t.Fatal(err)
	}
	return lt
}

func logicalTestLinearTransformation(t *testing.T, params ckks.Parameters, level int, diagonals map[int][]*big.Int, dimensions ring.Dimensions) lintrans.LinearTransformation {
	t.Helper()
	lt := lintrans.LinearTransformation{
		MetaData: &rlwe.MetaData{PlaintextMetaData: rlwe.PlaintextMetaData{
			Scale: rlwe.NewScale(1 << 20), IsBatched: true, LogDimensions: dimensions,
		}, CiphertextMetaData: rlwe.CiphertextMetaData{IsNTT: true, IsMontgomery: true}},
		LevelQ: level, N1: 0, Vec: make(map[int]ringqp.Poly, len(diagonals)),
	}
	for diagonal, values := range diagonals {
		lt.Vec[diagonal] = ringqp.Poly{Q: logicalQPolynomialFromIntegers(params, level, values)}
	}
	return lt
}

func decodeFastStoragePlaintext(t *testing.T, pt *FastStoragePlaintext) []*big.Int {
	t.Helper()
	coefficients := make([][]uint64, 3)
	for row := 0; row < 3; row++ {
		coefficients[row] = make([]uint64, pt.params.N())
		subring, _ := pt.basis.subring(row)
		subring.INTT(pt.value.Coeffs[row], coefficients[row])
	}
	values := make([]*big.Int, pt.params.N())
	for coefficient := range values {
		var residues [3]uint64
		for row := 0; row < 3; row++ {
			residues[row] = coefficients[row][coefficient]
		}
		value, err := pt.basis.decode(residues, 3)
		require.NoError(t, err)
		values[coefficient] = value
	}
	return values
}

func independentNegacyclicProduct(left, right []*big.Int, n int) []*big.Int {
	output := make([]*big.Int, n)
	for i := range output {
		output[i] = new(big.Int)
	}
	for i, x := range left {
		for j, y := range right {
			if x.Sign() == 0 || y.Sign() == 0 {
				continue
			}
			product := new(big.Int).Mul(x, y)
			index := i + j
			if index >= n {
				index -= n
				product.Neg(product)
			}
			output[index].Add(output[index], product)
		}
	}
	return output
}

func independentAutomorphism(input []*big.Int, n int, galEl uint64) []*big.Int {
	output := make([]*big.Int, n)
	for i := range output {
		output[i] = new(big.Int)
	}
	order := uint64(2 * n)
	for coefficient, value := range input {
		if value.Sign() == 0 {
			continue
		}
		exponent := (uint64(coefficient) * galEl) % order
		sign := int64(1)
		if exponent >= uint64(n) {
			exponent -= uint64(n)
			sign = -1
		}
		term := new(big.Int).Set(value)
		if sign < 0 {
			term.Neg(term)
		}
		output[exponent].Add(output[exponent], term)
	}
	return output
}
