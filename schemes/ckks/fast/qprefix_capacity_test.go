package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestObserveQPrefixCapacityReportsExactCenteredComponentBounds(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 40, 40, 45},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewEvaluator(params)
	const level, rows = 3, 4
	ct := NewCiphertext(params, 1, level)
	ct.IsNTT = false
	ct.Scale = rlwe.NewScale(1 << 40)
	componentMaxima := []int64{5, 7}
	for component := range ct.Value {
		for row := 0; row < rows; row++ {
			q := new(big.Int).SetUint64(params.Q()[row])
			for i := 0; i < params.N(); i++ {
				value := int64((i % int(2*componentMaxima[component]+1)) - int(componentMaxima[component]))
				ct.Value[component].Coeffs[row][i] = new(big.Int).Mod(big.NewInt(value), q).Uint64()
			}
		}
	}

	var got []QPrefixCapacitySnapshot
	eval.SetQPrefixCapacityObserver(func(snapshot QPrefixCapacitySnapshot) error {
		got = append(got, snapshot)
		return nil
	})
	require.NoError(t, eval.ObserveQPrefixCapacity("unit-checkpoint", ct, rows))
	require.Len(t, got, 1)
	prefixProduct := big.NewInt(1)
	for row := 0; row < rows; row++ {
		prefixProduct.Mul(prefixProduct, new(big.Int).SetUint64(params.Q()[row]))
	}
	require.Equal(t, QPrefixCapacitySnapshot{
		Name:          "unit-checkpoint",
		Level:         level,
		Rows:          rows,
		Scale:         ct.Scale.Value.Text('g', -1),
		Degree:        1,
		MaxAbs:        []string{"5", "7"},
		PrefixProduct: prefixProduct.String(),
		StrictFit:     true,
	}, got[0])
}
