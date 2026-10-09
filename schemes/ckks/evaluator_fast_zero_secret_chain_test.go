package ckks

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func requireFastChainOracle(t *testing.T, params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, ct *rlwe.Ciphertext, want []complex128, tolerance float64, checkpoint string) {
	t.Helper()
	got := decodeFastAddSubTestValues(t, params, encoder, decryptor, ct)
	var sumSquares float64
	var maxError float64
	for i := range want {
		error := cmplx.Abs(got[i] - want[i])
		sumSquares += error * error
		maxError = max(maxError, error)
	}
	rmse := math.Sqrt(sumSquares / float64(len(want)))
	require.LessOrEqual(t, rmse, tolerance, "%s RMSE", checkpoint)
	require.LessOrEqual(t, maxError, tolerance, "%s max complex error", checkpoint)
}

func TestPublicFastZeroSecretCompactAddMulRescaleRotateChain(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	decryptor := rlwe.NewDecryptor(params, sk)
	evaluator := NewEvaluator(params, nil)
	values0, values1 := fastAddSubTestValues(1), fastAddSubTestValues(5)

	for _, level := range []int{1, 3} {
		t.Run(fmt.Sprintf("level_%d", level), func(t *testing.T) {
			ct0 := encryptFastAddSubTestValues(t, params, encoder, encryptor, values0, level)
			ct1 := encryptFastAddSubTestValues(t, params, encoder, encryptor, values1, level)
			wantAdd := make([]complex128, len(values0))
			wantMul := make([]complex128, len(values0))
			for i := range wantAdd {
				wantAdd[i] = values0[i] + values1[i]
				wantMul[i] = wantAdd[i] * values1[i]
			}

			added, err := evaluator.AddNew(ct0, ct1)
			require.NoError(t, err)
			require.Equal(t, level, added.Level())
			require.True(t, ct0.Scale.Equal(added.Scale))
			requireFastChainOracle(t, params, encoder, decryptor, added, wantAdd, 1e-6, "Add")

			product, err := evaluator.MulRelinNew(added, ct1)
			require.NoError(t, err)
			require.Equal(t, level, product.Level())
			require.True(t, added.Scale.Mul(ct1.Scale).Equal(product.Scale))
			requireFastChainOracle(t, params, encoder, decryptor, product, wantMul, 1e-4, "MulRelin")

			rescaled := NewCiphertext(params, 1, level-1)
			require.NoError(t, evaluator.Rescale(product, rescaled))
			require.Equal(t, level-1, rescaled.Level())
			require.True(t, product.Scale.Div(rlwe.NewScale(params.Q()[level])).Equal(rescaled.Scale))
			requireFastChainOracle(t, params, encoder, decryptor, rescaled, wantMul, 1e-3, "Rescale")

			rotated, err := evaluator.RotateNew(rescaled, 1)
			require.NoError(t, err)
			require.Equal(t, level-1, rotated.Level())
			require.True(t, rescaled.Scale.Equal(rotated.Scale))
			requireFastChainOracle(t, params, encoder, decryptor, rotated, fastRotateTestOracle(wantMul, 1), 1e-3, "Rotate")
		})
	}
}
