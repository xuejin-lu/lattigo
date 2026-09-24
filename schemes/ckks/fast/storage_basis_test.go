package fast

import (
	"math/big"
	"math/bits"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func TestFastStorageBasisPrimesAndNTTRoundTrip(t *testing.T) {
	primes := fastStoragePrimes
	for i, q := range primes {
		require.True(t, ring.IsPrime(q))
		require.Equal(t, 60, bits.Len64(q))
		require.Less(t, q, uint64(1)<<60)
		require.Equal(t, uint64(1), q%(1<<18))
		for j := 0; j < i; j++ {
			require.NotEqual(t, primes[j], q)
		}
	}

	for _, logN := range []int{12, 13, 16} {
		basis, err := newFastStorageBasis(logN)
		require.NoError(t, err)
		wantProduct := big.NewInt(1)
		for _, q := range primes {
			wantProduct.Mul(wantProduct, new(big.Int).SetUint64(q))
		}
		gotProduct, err := basis.product(3)
		require.NoError(t, err)
		require.Equal(t, wantProduct, gotProduct)
		for i, q := range primes {
			subring, err := basis.subring(i)
			require.NoError(t, err)
			require.Equal(t, q, subring.Modulus)
			input := make([]uint64, subring.N)
			transformed := make([]uint64, subring.N)
			output := make([]uint64, subring.N)
			for j := range input {
				input[j] = uint64((j*j+17*j+31)%100003) % q
			}
			subring.NTT(input, transformed)
			subring.INTT(transformed, output)
			require.Equal(t, input, output, "LogN=%d primeIndex=%d", logN, i)
		}
	}
}

func TestFastStorageBasisFixedWidthSignedRoundTrip(t *testing.T) {
	basis, err := newFastStorageBasis(13)
	require.NoError(t, err)
	rng := rand.New(rand.NewSource(20260924))

	for width := 1; width <= 3; width++ {
		product, err := basis.product(width)
		require.NoError(t, err)
		capacity, err := basis.centeredCapacity(width)
		require.NoError(t, err)
		values := []*big.Int{
			big.NewInt(0), big.NewInt(1), big.NewInt(-1),
			new(big.Int).Set(capacity), new(big.Int).Neg(new(big.Int).Set(capacity)),
			new(big.Int).Lsh(big.NewInt(1), uint(width*59-1)),
		}
		values[len(values)-1].Mod(values[len(values)-1], capacity)
		for i := 0; i < 100; i++ {
			value := new(big.Int).Rand(rng, capacity)
			if i&1 != 0 {
				value.Neg(value)
			}
			values = append(values, value)
		}
		for _, value := range values {
			unique, err := basis.hasUniqueCenteredRepresentation(value, width)
			require.NoError(t, err)
			require.True(t, unique)
			encoded, err := basis.encode(value, width)
			require.NoError(t, err)
			for i := 0; i < width; i++ {
				wantResidue := new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(fastStoragePrimes[i]))
				require.Equal(t, wantResidue.Uint64(), encoded[i], "width=%d primeIndex=%d value=%s", width, i, value)
			}
			decoded, err := basis.decode(encoded, width)
			require.NoError(t, err)
			require.Equal(t, value, decoded, "width=%d value=%s product=%s", width, value, product)
			fixed, err := basis.decodeFixed(encoded, width)
			require.NoError(t, err)
			require.Equal(t, value, signedStorageBig(fixed))
		}
	}
}

func TestFastStorageBasisStrictCenteredCapacity(t *testing.T) {
	basis, err := newFastStorageBasis(13)
	require.NoError(t, err)
	for width := 1; width <= 3; width++ {
		product, err := basis.product(width)
		require.NoError(t, err)
		capacity, err := basis.centeredCapacity(width)
		require.NoError(t, err)
		require.Equal(t, new(big.Int).Rsh(new(big.Int).Sub(new(big.Int).Set(product), big.NewInt(1)), 1), capacity)
		ok, err := basis.hasUniqueCenteredRepresentation(capacity, width)
		require.NoError(t, err)
		require.True(t, ok)
		next := new(big.Int).Add(new(big.Int).Set(capacity), big.NewInt(1))
		ok, err = basis.hasUniqueCenteredRepresentation(next, width)
		require.NoError(t, err)
		require.False(t, ok)
		negNext := new(big.Int).Neg(new(big.Int).Set(next))
		positiveResidues, err := basis.encode(next, width)
		require.NoError(t, err)
		negativeResidues, err := basis.encode(negNext, width)
		require.NoError(t, err)
		decodedPositive, err := basis.decode(positiveResidues, width)
		require.NoError(t, err)
		decodedNegative, err := basis.decode(negativeResidues, width)
		require.NoError(t, err)
		require.Equal(t, new(big.Int).Neg(new(big.Int).Set(capacity)), decodedPositive)
		require.Equal(t, capacity, decodedNegative)
	}
}

func TestFastStorageBasisWidthValidationAndLogicalQIndependence(t *testing.T) {
	basis12, err := newFastStorageBasis(12)
	require.NoError(t, err)
	basis13, err := newFastStorageBasis(13)
	require.NoError(t, err)
	basis16, err := newFastStorageBasis(16)
	require.NoError(t, err)
	for width := 1; width <= 3; width++ {
		p12, err := basis12.product(width)
		require.NoError(t, err)
		p13, err := basis13.product(width)
		require.NoError(t, err)
		p16, err := basis16.product(width)
		require.NoError(t, err)
		require.Equal(t, p12, p13)
		require.Equal(t, p13, p16)
	}
	for _, width := range []int{0, 4, -1} {
		_, err := basis13.product(width)
		require.Error(t, err)
		_, err = basis13.decode([3]uint64{}, width)
		require.Error(t, err)
	}
	_, err = basis13.subring(3)
	require.Error(t, err)
	_, err = basis13.decode([3]uint64{fastStoragePrimes[0]}, 1)
	require.Error(t, err)
	tooWide := new(big.Int).Lsh(big.NewInt(1), 192)
	_, err = basis13.encode(tooWide, 3)
	require.Error(t, err)

	logicalGenerator := ring.NewNTTFriendlyPrimesGenerator(45, 1<<14)
	logicalQ, err := logicalGenerator.NextAlternatingPrimes(4)
	require.NoError(t, err)
	_, err = ring.NewRing(1<<13, logicalQ)
	require.NoError(t, err)
	for _, q := range logicalQ {
		for _, storagePrime := range fastStoragePrimes {
			require.NotEqual(t, storagePrime, q, "logical Q must remain separate from private storage basis")
		}
	}
}

func signedStorageBig(value storageInteger) *big.Int {
	result := storageToBig(value.magnitude)
	if value.negative {
		result.Neg(result)
	}
	return result
}
