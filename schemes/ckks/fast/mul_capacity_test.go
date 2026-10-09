package fast

import (
	"math"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestPublicFastMulRelinProductPassesExactQPrefixCapacityAtLevelsOneAndThree(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            8,
		LogQ:            []int{55, 39, 40, 39},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	encoder := ckks.NewEncoder(params)
	sk := ckks.NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	publicEval := ckks.NewEvaluator(params, nil)
	capacityEval := NewEvaluator(params)
	valuesA := []complex128{0.125 + 0.0625i, -0.25 + 0.125i, 0.0625 - 0.1875i, 0.3125, -0.125 - 0.0625i, 0.1875 + 0.25i, -0.0625, 0.25 - 0.125i}
	valuesB := []complex128{0.25 - 0.0625i, 0.125 + 0.1875i, -0.25 + 0.0625i, 0.0625 + 0.125i, 0.1875, -0.125 + 0.0625i, 0.3125 - 0.1875i, -0.0625 - 0.125i}

	for _, level := range []int{1, 3} {
		newInput := func(values []complex128) *rlwe.Ciphertext {
			pt := ckks.NewPlaintext(params, level)
			pt.Scale = rlwe.NewScale(math.Exp2(16))
			pt.IsNTT = true
			pt.LogDimensions = ring.Dimensions{Cols: params.LogMaxSlots()}
			require.NoError(t, encoder.Encode(values, pt))
			ct, err := encryptor.EncryptNew(pt)
			require.NoError(t, err)
			return ct
		}
		a, b := newInput(valuesA), newInput(valuesB)
		product, err := publicEval.MulRelinNew(a, b)
		require.NoError(t, err)
		var snapshots []QPrefixCapacitySnapshot
		capacityEval.SetQPrefixCapacityObserver(func(snapshot QPrefixCapacitySnapshot) error {
			snapshots = append(snapshots, snapshot)
			return nil
		})
		rows := level + 1
		require.NoError(t, capacityEval.ObserveQPrefixCapacity("public_mulrelin_output", product, rows))
		require.Len(t, snapshots, 1)
		require.True(t, snapshots[0].StrictFit, "Level %d component bounds must satisfy strict 2B < S_Q: %+v", level, snapshots[0])
		require.Equal(t, rows, snapshots[0].Rows)
		require.Equal(t, level, snapshots[0].Level)
	}
}
