package ckks

import (
	"fmt"
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
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
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	decryptor := rlwe.NewDecryptor(params, sk)

	for name, key := range map[string]rlwe.EncryptionKey{"SK": sk, "PK": pk} {
		t.Run(name, func(t *testing.T) {
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

						ct, err := rlwe.NewEncryptor(params, key).EncryptNew(pt)
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
		})
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

	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	for name, key := range map[string]rlwe.EncryptionKey{"SK": sk, "PK": pk} {
		t.Run(name, func(t *testing.T) {
			ct, err := rlwe.NewEncryptor(params, key).EncryptNew(pt)
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
		})
	}
}

// Exercise the same ordinary key lifecycle used by public-key applications.
func TestFastCKKSZeroSecretPublicKeyRotateAndHoisted(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	gk := kgen.GenGaloisKeyNew(params.GaloisElementForRotation(1), sk)
	eval := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
	encoder := NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	values := fastRotateTestValues(fastRotateTestLogSlots)
	input := fastRotateTestEncrypt(t, params, encoder, rlwe.NewEncryptor(params, pk), values, fastRotateTestLogSlots, 3, true)
	want := fastRotateTestOracle(values, 1)

	rotated, err := eval.RotateNew(input, 1)
	require.NoError(t, err)
	requireFastRotateOracle(t, params, encoder, decryptor, input, rotated, want)

	level := input.Level()
	decomposition := make([]ringqp.Poly, params.BaseRNSDecompositionVectorSize(level, gk.LevelP()))
	for i := range decomposition {
		decomposition[i] = params.RingQP().AtLevel(level, gk.LevelP()).NewPoly()
	}
	eval.DecomposeNTT(level, gk.LevelP(), gk.LevelP()+1, input.Value[1], input.IsNTT, decomposition)
	qp := rlwe.NewElementExtended(params, 1, level, gk.LevelP())
	require.NoError(t, eval.AutomorphismHoistedLazy(level, input, decomposition, gk.GaloisElement, qp))
	got := NewCiphertext(params, 1, level)
	*got.MetaData = *input.MetaData
	for i := range got.Value {
		eval.ModDownQPtoQNTT(level, gk.LevelP(), qp.Value[i].Q, qp.Value[i].P, got.Value[i])
	}
	requireFastRotateOracle(t, params, encoder, decryptor, input, got, want)
	require.Positive(t, eval.fastQPCompat.counterSnapshot().keySwitchElisions)
}

func TestFastCKKSZeroSecretEncryptContracts(t *testing.T) {
	params := fastZeroSecretEncryptionTestParameters(t)
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	pt := NewPlaintext(params, 1)
	pt.Value.Coeffs[0][0] = 7
	for name, key := range map[string]rlwe.EncryptionKey{"SK": sk, "PK": pk} {
		t.Run(name, func(t *testing.T) {
			enc := rlwe.NewEncryptor(params, key)
			out := NewCiphertext(params, 1, 0)
			out.Value[1].Coeffs[0][0] = 99
			require.NoError(t, enc.Encrypt(pt, out))
			require.Equal(t, 0, out.Level())
			require.Equal(t, pt.Value.Coeffs[0], out.Value[0].Coeffs[0])
			require.True(t, fastQPolyIsZero(out.Value[1], 1))
			require.True(t, pt.MetaData.Equal(out.MetaData))

			var nilCT *rlwe.Ciphertext
			require.ErrorContains(t, enc.Encrypt(pt, nilCT), "ciphertext cannot be nil")
			require.ErrorContains(t, enc.Encrypt(pt, NewCiphertext(params, 2, 1)), "degree-one storage")
			bad := NewCiphertext(params, 1, 1)
			bad.MetaData = nil
			require.ErrorContains(t, enc.Encrypt(pt, bad), "metadata must be non-nil")
			badPT := NewPlaintext(params, 1)
			badPT.MetaData = nil
			require.ErrorContains(t, enc.Encrypt(badPT, NewCiphertext(params, 1, 1)), "metadata must be non-nil")
			require.Error(t, enc.Encrypt(pt, &rlwe.Element[ring.Poly]{}))

			// Nil plaintext and EncryptZero retain native randomized encryption.
			for _, encryptZero := range []func(*rlwe.Ciphertext) error{
				func(ct *rlwe.Ciphertext) error { return enc.EncryptZero(ct) },
				func(ct *rlwe.Ciphertext) error { return enc.Encrypt(nil, ct) },
			} {
				zero := NewCiphertext(params, 1, 1)
				require.NoError(t, encryptZero(zero))
				require.False(t, fastQPolyIsZero(zero.Value[1], 2))
			}
		})
	}
	// WithKey must retain the capability and the existing public-key checks.
	ct, err := rlwe.NewEncryptor(params, sk).WithKey(pk).EncryptNew(pt)
	require.NoError(t, err)
	require.True(t, fastQPolyIsZero(ct.Value[1], 2))
	require.Error(t, rlwe.NewEncryptor(params, nil).Encrypt(pt, NewCiphertext(params, 1, 1)))
	invalid := pk.CopyNew()
	invalid.Layout = rlwe.KeyLayoutFast
	require.Panics(t, func() { rlwe.NewEncryptor(params, invalid) })
}
