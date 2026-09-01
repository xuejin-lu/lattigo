package fast

import (
	"bytes"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// Helper function to safely check if a ringqp.Poly is completely zero
func isZeroPoly(p ringqp.Poly) bool {
	// 直接走訪 Q 的係數（如果 p.Q.Coeffs 是空的，迴圈根本不會執行，非常安全）
	for _, coeffs := range p.Q.Coeffs {
		for _, c := range coeffs {
			if c != 0 {
				return false
			}
		}
	}

	// 直接走訪 P 的係數
	for _, coeffs := range p.P.Coeffs {
		for _, c := range coeffs {
			if c != 0 {
				return false
			}
		}
	}

	return true
}

func TestFastCKKSKeyGenerationPhase1B(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ExampleParameters128BitLogN14LogQP438)
	require.NoError(t, err)

	stdKgen := rlwe.NewKeyGenerator(params)
	fastKgen := NewKeyGenerator(params)

	sk := fastKgen.GenSecretKeyNew()

	t.Run("Fast PublicKey Behavior", func(t *testing.T) {
		pk := fastKgen.GenPublicKeyNew(sk)
		require.Equal(t, rlwe.KeyLayoutFast, pk.Layout)
		require.True(t, isZeroPoly(pk.Value[1]))
		require.False(t, isZeroPoly(pk.Value[0]))
		for _, limb := range pk.Value[0].Q.Coeffs {
			require.NotEmpty(t, limb)
		}
		for _, limb := range pk.Value[0].P.Coeffs {
			require.NotEmpty(t, limb)
		}
	})

	t.Run("Fast EvaluationKey material across all dimensions", func(t *testing.T) {
		evk, err := fastKgen.GenEvaluationKeyNew(sk, sk)
		require.NoError(t, err)
		require.Equal(t, rlwe.KeyLayoutFast, evk.Layout)
		require.Equal(t, rlwe.KeyLayoutFast, evk.GadgetCiphertext.Layout)
		for i := range evk.Value {
			for j := range evk.Value[i] {
				require.Len(t, evk.Value[i][j], 2)
				require.True(t, isZeroPoly(evk.Value[i][j][1]))
				require.False(t, isZeroPoly(evk.Value[i][j][0]))
				// The Fast construction calls standard EncryptZero with zero
				// secrets, hence first is exactly e and cannot contain a key term.
				for _, limb := range evk.Value[i][j][0].Q.Coeffs {
					require.NotEmpty(t, limb)
				}
				for _, limb := range evk.Value[i][j][0].P.Coeffs {
					require.NotEmpty(t, limb)
				}
			}
		}
	})

	t.Run("Fast Compressed Key Rejection", func(t *testing.T) {
		_, err := fastKgen.GenEvaluationKeyNew(sk, sk, rlwe.EvaluationKeyParameters{Compressed: true})
		require.Error(t, err)
	})

	t.Run("Fast Key EvaluationKey Expand Error", func(t *testing.T) {
		evk, err := fastKgen.GenEvaluationKeyNew(sk, sk)
		require.NoError(t, err)
		err = evk.Expand(params, nil)
		require.Error(t, err)
		require.Contains(t, err.Error(), "cannot be expanded")
	})

	t.Run("Serialization Rejection", func(t *testing.T) {
		pk := fastKgen.GenPublicKeyNew(sk)
		evk, err := fastKgen.GenEvaluationKeyNew(sk, sk)
		require.NoError(t, err)

		_, err = pk.MarshalBinary()
		require.Error(t, err)

		_, err = evk.MarshalBinary()
		require.Error(t, err)

		buf := new(bytes.Buffer)
		_, err = pk.WriteTo(buf)
		require.Error(t, err)
	})

	t.Run("Standard Consumer Rejection", func(t *testing.T) {
		pk := fastKgen.GenPublicKeyNew(sk)
		enc := rlwe.NewEncryptor(params, nil)
		require.Panics(t, func() { enc.WithKey(pk) })

		evk, err := fastKgen.GenEvaluationKeyNew(sk, sk)
		require.NoError(t, err)
		eval := rlwe.NewEvaluator(params, nil)
		err = eval.ApplyEvaluationKey(nil, evk, nil)
		require.Error(t, err)
		require.Contains(t, err.Error(), "Fast evaluation key")
	})

	t.Run("Standard Compatibility Regression", func(t *testing.T) {
		pk := stdKgen.GenPublicKeyNew(sk)
		require.Equal(t, rlwe.KeyLayoutStandard, pk.Layout)
		require.False(t, isZeroPoly(pk.Value[1])) // Normal 'a' is preserved

		_, err := pk.MarshalBinary()
		require.NoError(t, err) // Standard keys serialize flawlessly
	})
}
