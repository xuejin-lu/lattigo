package fast

import (
	"errors"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestQPrefixWidth(t *testing.T) {
	tests := []struct {
		level int
		width int
	}{
		{level: 0, width: 1},
		{level: 1, width: 2},
		{level: 2, width: 3},
		{level: 3, width: 4},
		{level: 4, width: 4},
		{level: 16, width: 4},
	}
	for _, test := range tests {
		got, err := QPrefixWidth(test.level)
		require.NoError(t, err)
		require.Equal(t, test.width, got, "level=%d", test.level)
	}
	_, err := QPrefixWidth(-1)
	require.Error(t, err)
}

func TestQPrefixProductUsesActualQValues(t *testing.T) {
	q := []uint64{101, 103, 107, 109, 113}
	tests := []struct {
		level int
		want  string
	}{
		{level: 0, want: "101"},
		{level: 1, want: "10403"},
		{level: 2, want: "1113121"},
		{level: 3, want: "121330189"},
		{level: 4, want: "121330189"},
	}
	for _, test := range tests {
		got, err := QPrefixProduct(q, test.level)
		require.NoError(t, err)
		require.Equal(t, test.want, got.String(), "level=%d", test.level)
	}
	_, err := QPrefixProduct(q[:3], 3)
	require.Error(t, err)
	_, err = QPrefixProduct([]uint64{101, 0, 107, 109}, 3)
	require.Error(t, err)
}

func TestQPrefixCapacityStrictBoundary(t *testing.T) {
	product := big.NewInt(10)
	require.True(t, QPrefixCapacitySatisfied(big.NewInt(4), product), "8 < 10 must fit")
	require.False(t, QPrefixCapacitySatisfied(big.NewInt(5), product), "the exact 2B == S_Q boundary must fail")
	require.False(t, QPrefixCapacitySatisfied(big.NewInt(6), product))
	require.False(t, QPrefixCapacitySatisfied(big.NewInt(-1), product))
	require.False(t, QPrefixCapacitySatisfied(big.NewInt(0), big.NewInt(0)))
}

func TestCheckQPrefixCapacityReportsComponent(t *testing.T) {
	q := []uint64{101, 103, 107, 109}
	product, err := QPrefixProduct(q, 2)
	require.NoError(t, err)

	require.NoError(t, CheckQPrefixCapacity(2, q, []*big.Int{
		big.NewInt(10),
		new(big.Int).Sub(new(big.Int).Rsh(new(big.Int).Set(product), 1), big.NewInt(1)),
	}))

	bound := new(big.Int).Rsh(new(big.Int).Set(product), 1)
	bound.Add(bound, big.NewInt(1))
	err = CheckQPrefixCapacity(2, q, []*big.Int{big.NewInt(10), bound})
	require.Error(t, err)
	var capacityErr *QPrefixCapacityError
	require.ErrorAs(t, err, &capacityErr)
	require.Equal(t, 2, capacityErr.Level)
	require.Equal(t, 1, capacityErr.Component)
	require.Equal(t, bound, capacityErr.Bound)
	require.Equal(t, product, capacityErr.PrefixProduct)
	require.Contains(t, err.Error(), "strict 2B < S_Q failed")

	before := new(big.Int).Set(bound)
	_ = CheckQPrefixCapacity(2, q, []*big.Int{bound})
	require.Equal(t, before, bound, "capacity preflight must not mutate caller bounds")
}

func TestCommitQPrefixBoundsIsTransactional(t *testing.T) {
	q := []uint64{101, 103, 107, 109}
	dst := []*big.Int{big.NewInt(3), big.NewInt(7)}
	original0 := new(big.Int).Set(dst[0])
	original1 := new(big.Int).Set(dst[1])
	originalFirst := dst[0]
	originalSecond := dst[1]

	product, err := QPrefixProduct(q, 1)
	require.NoError(t, err)
	tooLarge := new(big.Int).Set(product)
	err = CommitQPrefixBounds(&dst, 1, q, []*big.Int{big.NewInt(5), tooLarge})
	require.Error(t, err)
	var capacityErr *QPrefixCapacityError
	require.ErrorAs(t, err, &capacityErr)
	require.Same(t, originalFirst, dst[0])
	require.Same(t, originalSecond, dst[1])
	require.Equal(t, original0, dst[0])
	require.Equal(t, original1, dst[1])

	candidates := []*big.Int{big.NewInt(11), big.NewInt(13)}
	require.NoError(t, CommitQPrefixBounds(&dst, 1, q, candidates))
	require.Equal(t, []string{"11", "13"}, []string{dst[0].String(), dst[1].String()})
	dst[0].SetInt64(99)
	require.Equal(t, "11", candidates[0].String(), "committed values must not alias candidates")

	priorFirst := dst[0]
	priorSecond := dst[1]
	prior0 := new(big.Int).Set(dst[0])
	prior1 := new(big.Int).Set(dst[1])
	err = CommitQPrefixBounds(&dst, 0, q, []*big.Int{nil})
	require.Error(t, err)
	require.Same(t, priorFirst, dst[0])
	require.Same(t, priorSecond, dst[1])
	require.Equal(t, prior0, dst[0])
	require.Equal(t, prior1, dst[1], "invalid metadata must not partially replace the destination")
}

func TestQPrefixCapacityErrorSupportsErrorsAs(t *testing.T) {
	q := []uint64{101}
	err := CheckQPrefixCapacity(0, q, []*big.Int{big.NewInt(51)})
	var capacityErr *QPrefixCapacityError
	require.True(t, errors.As(err, &capacityErr))
}
