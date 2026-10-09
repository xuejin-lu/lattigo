package ckks

import (
	"fmt"
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func fastZeroSecretEncryptionTestParameters(t *testing.T) Parameters {
	t.Helper()
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN:            8,
		LogQ:            []int{55, 45},
		LogDefaultScale: 40,
	})
	require.NoError(t, err)
	return params
}

func TestFastCKKSZeroSecretEncryptNewPreservesPlaintext(t *testing.T) {
	params := fastZeroSecretEncryptionTestParameters(t)
	encoder := NewEncoder(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	decryptor := rlwe.NewDecryptor(params, sk)

	for _, level := range []int{0, 1} {
		for _, isNTT := range []bool{false, true} {
			t.Run(fmt.Sprintf("level=%d/ntt=%t", level, isNTT), func(t *testing.T) {
				logSlots := params.LogMaxSlots()
				if isNTT {
					logSlots = 4
				}
				values := make([]complex128, 1<<logSlots)
				for i := range values {
					values[i] = complex(float64(i%17-8)/32, float64((i%5)-2)/64)
				}
				pt := NewPlaintext(params, level)
				pt.IsNTT = isNTT
				pt.LogDimensions = ring.Dimensions{Cols: logSlots}
				require.NoError(t, encoder.Encode(values, pt))
				originalPoly := pt.Value.CopyNew()

				ct, err := rlwe.NewEncryptor(params, sk).EncryptNew(pt)
				require.NoError(t, err)
				require.Equal(t, 1, ct.Degree())
				require.Equal(t, level, ct.Level())
				require.True(t, pt.MetaData.Equal(ct.MetaData))
				require.NotSame(t, pt.MetaData, ct.MetaData)
				for row := 0; row <= level; row++ {
					require.Equal(t, pt.Value.Coeffs[row], ct.Value[0].Coeffs[row])
					for _, coefficient := range ct.Value[1].Coeffs[row] {
						require.Zero(t, coefficient)
					}
				}

				decoded := make([]complex128, len(values))
				require.NoError(t, encoder.Decode(decryptor.DecryptNew(ct), decoded))
				for i := range values {
					require.LessOrEqual(t, math.Abs(real(decoded[i]-values[i])), 1e-7)
					require.LessOrEqual(t, math.Abs(imag(decoded[i]-values[i])), 1e-7)
				}

				ct.Value[0].Coeffs[0][0] ^= 1
				require.True(t, pt.Value.Equal(originalPoly), "mutating the output ciphertext must not mutate the input plaintext")
			})
		}
	}
}

func TestFastCKKSZeroSecretEncryptNewCopiesMetadataDeeply(t *testing.T) {
	params := fastZeroSecretEncryptionTestParameters(t)
	pt := NewPlaintext(params, 1)
	// Synthetic flags exercise metadata preservation only; this fixture is not
	// used for numerical decoding in these non-default metadata states.
	pt.IsNTT = false
	pt.IsMontgomery = true
	pt.IsBatched = false
	pt.IsBitReversed = true
	pt.LogDimensions = ring.Dimensions{Rows: 1, Cols: 3}
	pt.Scale = rlwe.NewScaleModT(37, 257)

	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	ct, err := rlwe.NewEncryptor(params, sk).EncryptNew(pt)
	require.NoError(t, err)
	require.True(t, pt.MetaData.Equal(ct.MetaData))
	require.NotSame(t, pt.MetaData, ct.MetaData)

	ptScale := pt.Scale.Float64()
	ptMod := new(big.Int).Set(pt.Scale.Mod)
	ct.Scale.Value.SetInt64(7)
	ct.Scale.Mod.SetInt64(11)
	ct.LogDimensions.Cols = 9
	require.Equal(t, ptScale, pt.Scale.Float64())
	require.Equal(t, ptMod, pt.Scale.Mod)
	require.Equal(t, 3, pt.LogDimensions.Cols)
}
