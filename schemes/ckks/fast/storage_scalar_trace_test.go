package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFastStorageMulIntegerSignedZeroLargeAndBounds(t *testing.T) {
	params := storageBoundaryParameters(t)
	values := [][]*big.Int{{big.NewInt(7), big.NewInt(-11), big.NewInt(19)}, {big.NewInt(-23), big.NewInt(31)}}
	large := new(big.Int).Add(new(big.Int).Lsh(big.NewInt(1), 70), big.NewInt(17))
	for _, ntt := range []bool{false, true} {
		for _, scalar := range []*big.Int{big.NewInt(3), big.NewInt(-5), large, new(big.Int).Neg(new(big.Int).Set(large)), new(big.Int)} {
			input := newFastStorageCiphertextFromBig(t, params, values, 2, 3, ntt, nil)
			input.metadata.Scale = input.Scale().Mul(rlwe.NewScale(13))
			input.metadata.IsBatched = false
			input.metadata.IsBitReversed = true
			input.metadata.LogDimensions = ring.Dimensions{Rows: 2, Cols: 3}
			before := snapshotFastStorageCiphertext(input)
			want := make([][]*big.Int, len(values))
			for component := range values {
				want[component] = make([]*big.Int, params.N())
				for coefficient, value := range values[component] {
					want[component][coefficient] = new(big.Int).Mul(value, scalar)
				}
				for coefficient := len(values[component]); coefficient < params.N(); coefficient++ {
					want[component][coefficient] = new(big.Int)
				}
			}

			output, err := FastStorageMulInteger(input, scalar)
			require.NoError(t, err)
			requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
			require.Equal(t, 3, output.StorageWidth())
			require.Equal(t, input.LogicalLevel(), output.LogicalLevel())
			require.Equal(t, ntt, output.IsNTT())
			require.False(t, output.IsMontgomery())
			require.True(t, input.Scale().Equal(output.Scale()), "integer multiplication must leave Scale unchanged")
			require.Equal(t, input.MetaData().IsBatched, output.MetaData().IsBatched)
			require.Equal(t, input.MetaData().IsBitReversed, output.MetaData().IsBitReversed)
			require.Equal(t, input.MetaData().LogDimensions, output.MetaData().LogDimensions)
			for component, bound := range input.ComponentBounds() {
				wantBound := new(big.Int).Mul(bound, new(big.Int).Abs(new(big.Int).Set(scalar)))
				require.Zero(t, wantBound.Cmp(output.ComponentBounds()[component]))
			}
			require.Equal(t, before, snapshotFastStorageCiphertext(input), "scalar multiplication must not mutate its input")
		}
	}
}

func TestFastStorageMulIntegerCapacityFailureIsTransactional(t *testing.T) {
	params := storageBoundaryParameters(t)
	input := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{new(big.Int)}}, 1, 3, true, nil)
	input.componentBounds[0] = new(big.Int).Rsh(storageToBig(input.basis.products[3]), 2)
	require.NoError(t, validateFastStorageState(input))
	before := snapshotFastStorageCiphertext(input)
	output, err := FastStorageMulInteger(input, big.NewInt(3))
	require.Nil(t, output)
	require.ErrorContains(t, err, "capacity")
	require.Equal(t, before, snapshotFastStorageCiphertext(input))
}

func TestFastStorageTraceNormalizedMonomialTheoremAndDivisibility(t *testing.T) {
	for _, ringLogN := range []int{4, 5} {
		t.Run("LogN-"+big.NewInt(int64(ringLogN)).String(), func(t *testing.T) {
			params := storageBoundaryParametersLogN(t, ringLogN)
			basis, err := fastStorageBasisForLogN(ringLogN)
			require.NoError(t, err)
			for _, logSlots := range []int{ringLogN - 1, ringLogN - 2, 0} {
				for _, exponent := range []int{0, 1, params.N() / 2, params.N() - 1} {
					t.Run("logSlots-"+big.NewInt(int64(logSlots)).String()+"-monomial-"+big.NewInt(int64(exponent)).String(), func(t *testing.T) {
						values := [][]*big.Int{storageTestMonomial(params.N(), exponent, 1), storageTestMonomial(params.N(), exponent, -1)}
						input := newFastStorageTraceCiphertext(t, params, basis, values)
						input.metadata.LogDimensions = ring.Dimensions{Rows: 1, Cols: 4}
						before := snapshotFastStorageCiphertext(input)
						wantSum, gap := independentStorageTrace(params, values, logSlots)
						want := make([][]*big.Int, len(wantSum))
						for component := range wantSum {
							want[component] = make([]*big.Int, params.N())
							nonZero := 0
							for coefficient, value := range wantSum[component] {
								quotient, remainder := new(big.Int), new(big.Int)
								quotient.QuoRem(value, big.NewInt(int64(gap)), remainder)
								require.Zero(t, remainder.Sign(), "Trace sum must divide by gap=%d", gap)
								want[component][coefficient] = quotient
								if quotient.Sign() != 0 {
									nonZero++
									require.LessOrEqual(t, new(big.Int).Abs(new(big.Int).Set(quotient)).Cmp(big.NewInt(1)), 0, "normalized monomial coefficient must be 0 or +/-1")
								}
							}
							require.LessOrEqual(t, nonZero, 1, "a monomial Trace must vanish or map to one signed monomial")
						}

						output, err := FastStorageTraceNormalized(input, logSlots)
						require.NoError(t, err)
						requireBigIntMatricesEqual(t, want, decodeFastStoragePolynomials(t, output))
						require.Equal(t, 3, output.StorageWidth())
						require.Equal(t, input.LogicalLevel(), output.LogicalLevel())
						require.True(t, input.Scale().Equal(output.Scale()))
						require.True(t, output.IsNTT())
						require.False(t, output.IsMontgomery())
						require.Equal(t, input.MetaData().LogDimensions, output.MetaData().LogDimensions)
						require.Equal(t, input.ComponentBounds(), output.ComponentBounds())
						require.Equal(t, before, snapshotFastStorageCiphertext(input), "Trace must not mutate its input")
					})
				}
			}
		})
	}
}

func TestFastStorageTraceNormalizedCapacityFailureIsTransactional(t *testing.T) {
	params := storageBoundaryParametersLogN(t, 4)
	input := newFastStorageCiphertextFromBig(t, params, [][]*big.Int{{new(big.Int)}, {new(big.Int)}}, 1, 3, true, nil)
	gap := 1 << uint(params.LogN()-1-1)
	capacity := storageToBig(input.basis.products[3])
	bound := new(big.Int).Quo(capacity, new(big.Int).SetUint64(uint64(2*gap)))
	bound.Add(bound, big.NewInt(1))
	input.componentBounds[0].Set(bound)
	input.componentBounds[1].Set(bound)
	require.NoError(t, validateFastStorageState(input), "fixture must fit before Trace")
	before := snapshotFastStorageCiphertext(input)
	output, err := FastStorageTraceNormalized(input, 1)
	require.Nil(t, output)
	require.ErrorContains(t, err, "strict 2*g*B < S3")
	require.Equal(t, before, snapshotFastStorageCiphertext(input), "capacity rejection must not mutate the input")
}

func newFastStorageTraceCiphertext(t *testing.T, params ckks.Parameters, basis fastStorageBasis, polynomials [][]*big.Int) *FastCiphertext {
	t.Helper()
	ct := newFastCiphertextWithBasis(params, len(polynomials)-1, 2, 3, basis)
	ct.metadata.IsNTT = true
	ct.componentBounds = zeroComponentBounds(len(polynomials))
	for component := range polynomials {
		for coefficient, value := range polynomials[component] {
			integer := storageIntegerFromBig(value)
			encoded, err := encodeCenteredStorageValue(basis, integer, 3)
			require.NoError(t, err)
			for row := 0; row < 3; row++ {
				ct.value[component].Coeffs[row][coefficient] = encoded[row]
			}
			magnitude := new(big.Int).Abs(new(big.Int).Set(value))
			if magnitude.Cmp(ct.componentBounds[component]) > 0 {
				ct.componentBounds[component].Set(magnitude)
			}
		}
		for row := 0; row < 3; row++ {
			subring, err := basis.subring(row)
			require.NoError(t, err)
			subring.NTT(ct.value[component].Coeffs[row], ct.value[component].Coeffs[row])
		}
	}
	require.NoError(t, validateFastStorageState(ct))
	return ct
}

func storageTestMonomial(n, exponent int, coefficient int64) []*big.Int {
	polynomial := make([]*big.Int, n)
	for i := range polynomial {
		polynomial[i] = new(big.Int)
	}
	polynomial[exponent].SetInt64(coefficient)
	return polynomial
}

func independentStorageTrace(params ckks.Parameters, input [][]*big.Int, logN int) ([][]*big.Int, int) {
	current := cloneStorageBigIntMatrix(input)
	gap := 1 << uint(params.LogN()-logN-1)
	if logN == 0 {
		gap <<= 1
	}
	apply := func(galEl uint64) {
		nthRoot := uint64(2 * params.N())
		for component := range current {
			permuted := make([]*big.Int, params.N())
			for i := range permuted {
				permuted[i] = new(big.Int)
			}
			for exponent, coefficient := range current[component] {
				mapped := uint64(exponent) * galEl % nthRoot
				if mapped >= uint64(params.N()) {
					permuted[mapped-uint64(params.N())].Sub(permuted[mapped-uint64(params.N())], coefficient)
				} else {
					permuted[mapped].Add(permuted[mapped], coefficient)
				}
			}
			for i := range current[component] {
				current[component][i].Add(current[component][i], permuted[i])
			}
		}
	}
	for i := logN; i < params.LogN()-1; i++ {
		apply(params.GaloisElement(1 << uint(i)))
	}
	if logN == 0 {
		apply(uint64(2*params.N()) - 1)
	}
	return current, gap
}

func cloneStorageBigIntMatrix(values [][]*big.Int) [][]*big.Int {
	copy := make([][]*big.Int, len(values))
	for component := range values {
		copy[component] = make([]*big.Int, len(values[component]))
		for coefficient, value := range values[component] {
			copy[component][coefficient] = new(big.Int).Set(value)
		}
	}
	return copy
}
