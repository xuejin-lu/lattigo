package ckks

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

const (
	fastRotateTestLogSlots = 4
	fastRotateTestLogScale = 30
)

func fastRotateTestParameters(t *testing.T, ringType ring.Type) Parameters {
	t.Helper()
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN:            8,
		LogQ:            []int{55, 45, 40, 39, 38, 37, 36},
		LogP:            []int{55},
		LogDefaultScale: fastRotateTestLogScale,
		RingType:        ringType,
	})
	require.NoError(t, err)
	return params
}

func fastRotateTestValues(logSlots int) []complex128 {
	values := make([]complex128, 1<<logSlots)
	for i := range values {
		values[i] = complex(float64(i+1)/64, float64((i*7)%11-5)/128)
	}
	return values
}

func fastRotateTestOracle(values []complex128, k int) []complex128 {
	rotated := make([]complex128, len(values))
	shift := k % len(values)
	if shift < 0 {
		shift += len(values)
	}
	for i := range rotated {
		rotated[i] = values[(i+shift)%len(values)]
	}
	return rotated
}

func fastRotateTestEncrypt(t *testing.T, params Parameters, encoder *Encoder, encryptor *rlwe.Encryptor, values []complex128, logSlots, level int, isNTT bool) *rlwe.Ciphertext {
	t.Helper()
	pt := NewPlaintext(params, level)
	pt.Scale = rlwe.NewScale(math.Exp2(fastRotateTestLogScale))
	pt.IsNTT = isNTT
	pt.LogDimensions = ring.Dimensions{Cols: logSlots}
	require.NoError(t, encoder.Encode(values, pt))
	ct, err := encryptor.EncryptNew(pt)
	require.NoError(t, err)
	require.Equal(t, 1, ct.Degree())
	require.Equal(t, level, ct.Level())
	require.Equal(t, isNTT, ct.IsNTT)
	require.False(t, ct.IsMontgomery)
	for row := 0; row <= level; row++ {
		for _, coefficient := range ct.Value[1].Coeffs[row] {
			require.Zero(t, coefficient, "Fast EncryptNew must produce c1=0 at q%d", row)
		}
	}
	return ct
}

func requireFastRotateOracle(t *testing.T, params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, in, got *rlwe.Ciphertext, expected []complex128) {
	t.Helper()
	require.Equal(t, 1, got.Degree())
	require.Equal(t, in.Level(), got.Level())
	require.True(t, in.Scale.Equal(got.Scale))
	require.Equal(t, in.IsNTT, got.IsNTT)
	require.Equal(t, in.IsMontgomery, got.IsMontgomery)
	require.Equal(t, *in.MetaData, *got.MetaData)
	for component := range got.Value {
		require.Len(t, got.Value[component].Coeffs, got.Level()+1)
		for row := 0; row <= got.Level(); row++ {
			require.Len(t, got.Value[component].Coeffs[row], params.N())
		}
	}
	for row := 0; row <= got.Level(); row++ {
		for _, coefficient := range got.Value[1].Coeffs[row] {
			require.Zero(t, coefficient, "rotated zero-secret c1 must remain zero at q%d", row)
		}
	}
	decoded := make([]complex128, len(expected))
	require.NoError(t, encoder.Decode(decryptor.DecryptNew(got), decoded))
	for i := range expected {
		require.LessOrEqual(t, cmplx.Abs(decoded[i]-expected[i]), 1e-6, "slot %d", i)
	}
}

type recordingPublicAutomorphismCore struct {
	fastcore.AutomorphismCore
	rows      []int
	pairCount []int
}

func (core *recordingPublicAutomorphismCore) ApplyRows(ringQ *ring.Ring, galEl uint64, isNTT bool, rows int, pairs ...fastcore.PolynomialPair) error {
	core.rows = append(core.rows, rows)
	core.pairCount = append(core.pairCount, len(pairs))
	return core.AutomorphismCore.ApplyRows(ringQ, galEl, isNTT, rows, pairs...)
}

func TestPublicFastZeroSecretRotateUsesSharedCoreForAllActiveRows(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	decryptor := rlwe.NewDecryptor(params, sk)
	eval := NewEvaluator(params, nil)
	recorder := &recordingPublicAutomorphismCore{AutomorphismCore: eval.automorphismCore}
	eval.automorphismCore = recorder
	values := fastRotateTestValues(fastRotateTestLogSlots)

	for _, level := range []int{1, 5} {
		for _, isNTT := range []bool{true} {
			input := fastRotateTestEncrypt(t, params, encoder, encryptor, values, fastRotateTestLogSlots, level, isNTT)
			if level == 5 {
				// Dormant q4/q5 are deliberately malformed and nonzero: neither
				// the public adapter nor the shared core may inspect those rows.
				for component := range input.Value {
					for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
						for i := range input.Value[component].Coeffs[row] {
							input.Value[component].Coeffs[row][i] = ^uint64(0)
						}
					}
				}
			}
			for _, k := range []int{0, 1, -1} {
				t.Run(fmt.Sprintf("level=%d/ntt=%t/k=%d", level, isNTT, k), func(t *testing.T) {
					got, err := eval.RotateNew(input, k)
					require.NoError(t, err)
					if level <= 3 {
						requireFastRotateOracle(t, params, encoder, decryptor, input, got, fastRotateTestOracle(values, k))
					} else {
						require.Equal(t, 1, got.Degree())
						require.Equal(t, level, got.Level())
						for component := range got.Value {
							for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
								require.Len(t, got.Value[component].Coeffs[row], params.N())
							}
							for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
								require.Nil(t, got.Value[component].Coeffs[row])
							}
						}
					}
					require.Equal(t, min(level+1, fastcore.MaxQPrefixWidth), recorder.rows[len(recorder.rows)-1], "public P0 must pass only maintained Q-prefix rows")
					require.Equal(t, 2, recorder.pairCount[len(recorder.pairCount)-1], "the common core must process c0 and c1")

					out := NewCiphertext(params, 1, level)
					out.IsNTT, out.IsMontgomery = input.IsNTT, input.IsMontgomery
					require.NoError(t, eval.Rotate(input, k, out))
					if level <= 3 {
						requireFastRotateOracle(t, params, encoder, decryptor, input, out, fastRotateTestOracle(values, k))
					} else {
						for component := range out.Value {
							for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
								require.Nil(t, out.Value[component].Coeffs[row])
							}
						}
					}

					if k == 1 {
						alias := input.CopyNew()
						require.NoError(t, eval.Rotate(alias, k, alias))
						if level <= 3 {
							requireFastRotateOracle(t, params, encoder, decryptor, input, alias, fastRotateTestOracle(values, k))
						} else {
							for component := range alias.Value {
								for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
									require.Nil(t, alias.Value[component].Coeffs[row])
								}
							}
						}
					}
				})
			}
		}
	}
}

func TestPublicFastZeroSecretRotateCoefficientDomainAtFullSlots(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	values := fastRotateTestValues(params.LogMaxSlots())
	input := fastRotateTestEncrypt(t, params, encoder, encryptor, values, params.LogMaxSlots(), 5, false)
	eval := NewEvaluator(params, nil)
	for component := range input.Value {
		for row := fastcore.MaxQPrefixWidth; row <= input.Level(); row++ {
			for i := range input.Value[component].Coeffs[row] {
				input.Value[component].Coeffs[row][i] = ^uint64(0)
			}
		}
	}

	got, err := eval.RotateNew(input, 1)
	require.NoError(t, err)
	require.Equal(t, input.Level(), got.Level())
	require.Equal(t, 1, got.Degree())
	require.False(t, got.IsNTT)
	for component := range got.Value {
		for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
			require.Len(t, got.Value[component].Coeffs[row], params.N())
		}
		for row := fastcore.MaxQPrefixWidth; row <= got.Level(); row++ {
			require.Nil(t, got.Value[component].Coeffs[row])
		}
	}
}

type rotateTrackingEvaluationKeySet struct {
	rlwe.EvaluationKeySet
	galoisKeyGets int
}

func (keys *rotateTrackingEvaluationKeySet) GetGaloisKey(galEl uint64) (*rlwe.GaloisKey, error) {
	keys.galoisKeyGets++
	return keys.EvaluationKeySet.GetGaloisKey(galEl)
}

func TestPublicFastZeroSecretRotateDoesNotAccessAvailableGaloisKey(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	galEl := params.GaloisElementForRotation(1)
	galoisKey := NewKeyGenerator(params).GenGaloisKeyNew(galEl, sk)
	require.Equal(t, rlwe.KeyLayoutStandard, galoisKey.Layout)
	keys := &rotateTrackingEvaluationKeySet{EvaluationKeySet: rlwe.NewMemEvaluationKeySet(nil, galoisKey)}
	eval := NewEvaluator(params, keys)
	keys.galoisKeyGets = 0 // Ignore constructor key-list enumeration.

	encoder := NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	values := fastRotateTestValues(fastRotateTestLogSlots)
	input := fastRotateTestEncrypt(t, params, encoder, rlwe.NewEncryptor(params, sk), values, fastRotateTestLogSlots, 3, true)
	got, err := eval.RotateNew(input, 1)
	require.NoError(t, err)
	require.Zero(t, keys.galoisKeyGets)
	requireFastRotateOracle(t, params, encoder, decryptor, input, got, fastRotateTestOracle(values, 1))
}

func TestPublicFastZeroSecretRotateRejectsInvalidInputsTransactionally(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	galEl := params.GaloisElementForRotation(1)
	galoisKey := NewKeyGenerator(params).GenGaloisKeyNew(galEl, sk)
	keys := &rotateTrackingEvaluationKeySet{EvaluationKeySet: rlwe.NewMemEvaluationKeySet(nil, galoisKey)}
	require.Equal(t, rlwe.KeyLayoutStandard, galoisKey.Layout)
	eval := NewEvaluator(params, keys)
	keys.galoisKeyGets = 0 // Ignore constructor key-list enumeration.
	encryptor := rlwe.NewEncryptor(params, sk)
	values := fastRotateTestValues(fastRotateTestLogSlots)
	valid := fastRotateTestEncrypt(t, params, encoder, encryptor, values, fastRotateTestLogSlots, 5, true)

	t.Run("nonzero-c1", func(t *testing.T) {
		input := valid.CopyNew()
		input.Value[1].Coeffs[2][3] = 1
		inputBefore := input.CopyNew()
		out := NewCiphertext(params, 1, 5)
		out.IsNTT = true
		out.Scale = rlwe.NewScale(123)
		outBefore := out.CopyNew()
		require.Error(t, eval.Rotate(input, 1, out))
		require.True(t, input.Equal(inputBefore))
		require.True(t, out.Equal(outBefore))
		require.Zero(t, keys.galoisKeyGets, "rejected Fast input must not consult an available Standard-layout Galois key")
	})

	t.Run("compact-level-5", func(t *testing.T) {
		input := valid.CopyNew()
		for component := range input.Value {
			input.Value[component].Coeffs[4] = nil
			input.Value[component].Coeffs[5] = nil
		}
		out, err := eval.RotateNew(input, 1)
		require.NoError(t, err)
		require.Equal(t, 5, out.Level())
		for component := range out.Value {
			for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
				require.Len(t, out.Value[component].Coeffs[row], params.N())
			}
			for row := fastcore.MaxQPrefixWidth; row <= out.Level(); row++ {
				require.Nil(t, out.Value[component].Coeffs[row], "dormant q%d must remain unmaterialized", row)
			}
		}
	})

	t.Run("missing-active-input-row", func(t *testing.T) {
		input := valid.CopyNew()
		input.Value[0].Coeffs[3] = nil
		out := NewCiphertext(params, 1, 5)
		out.IsNTT = true
		outBefore := out.CopyNew()
		require.Error(t, eval.Rotate(input, 1, out))
		require.Nil(t, input.Value[0].Coeffs[3])
		for component := range input.Value {
			for row := range input.Value[component].Coeffs {
				if component == 0 && row == 3 {
					continue
				}
				require.Equal(t, valid.Value[component].Coeffs[row], input.Value[component].Coeffs[row])
			}
		}
		require.True(t, out.Equal(outBefore))
	})

	t.Run("noncanonical-c0-residue", func(t *testing.T) {
		input := valid.CopyNew()
		input.Value[0].Coeffs[2][3] = params.RingQ().SubRings[2].Modulus
		inputBefore := input.CopyNew()
		out := NewCiphertext(params, 1, 5)
		out.IsNTT = true
		outBefore := out.CopyNew()
		require.Error(t, eval.Rotate(input, 1, out))
		require.Equal(t, inputBefore.Value[0].Coeffs, input.Value[0].Coeffs)
		require.True(t, out.Equal(outBefore))
	})

	t.Run("missing-output-row", func(t *testing.T) {
		out := NewCiphertext(params, 1, 5)
		out.Value[0].Coeffs[3] = nil
		outBefore := out.CopyNew()
		require.Error(t, eval.Rotate(valid, 1, out))
		require.True(t, out.Equal(outBefore))
	})

	t.Run("mismatched-domain", func(t *testing.T) {
		out := NewCiphertext(params, 1, 5)
		out.IsNTT = false
		outBefore := out.CopyNew()
		require.Error(t, eval.Rotate(valid, 1, out))
		require.True(t, out.Equal(outBefore))
	})

	t.Run("invalid-montgomery-coefficient-state", func(t *testing.T) {
		input, out := valid.CopyNew(), NewCiphertext(params, 1, 5)
		input.IsNTT, input.IsMontgomery = false, true
		out.IsNTT, out.IsMontgomery = false, true
		inputBefore, outBefore := input.CopyNew(), out.CopyNew()
		require.Error(t, eval.Rotate(input, 1, out))
		require.True(t, input.Equal(inputBefore))
		require.True(t, out.Equal(outBefore))
	})

	t.Run("nil-metadata", func(t *testing.T) {
		input, out := valid.CopyNew(), NewCiphertext(params, 1, 5)
		out.IsNTT = true
		input.MetaData = nil
		outBefore := out.CopyNew()
		inputC0, inputC1 := input.Value[0].CopyNew(), input.Value[1].CopyNew()
		require.Error(t, eval.Rotate(input, 1, out))
		require.Nil(t, input.MetaData)
		require.Equal(t, inputC0.Coeffs, input.Value[0].Coeffs)
		require.Equal(t, inputC1.Coeffs, input.Value[1].Coeffs)
		require.True(t, out.Equal(outBefore))
	})

	t.Run("nil-output-metadata", func(t *testing.T) {
		out := NewCiphertext(params, 1, 5)
		out.IsNTT = true
		out.MetaData = nil
		outC0, outC1 := out.Value[0].CopyNew(), out.Value[1].CopyNew()
		require.Error(t, eval.Rotate(valid, 1, out))
		require.Nil(t, out.MetaData)
		require.Equal(t, outC0.Coeffs, out.Value[0].Coeffs)
		require.Equal(t, outC1.Coeffs, out.Value[1].Coeffs)
	})

	t.Run("degree-two", func(t *testing.T) {
		input, out := NewCiphertext(params, 2, 5), NewCiphertext(params, 1, 5)
		input.IsNTT, out.IsNTT = true, true
		inputBefore, outBefore := input.CopyNew(), out.CopyNew()
		require.Error(t, eval.Rotate(input, 1, out))
		require.True(t, input.Equal(inputBefore))
		require.True(t, out.Equal(outBefore))
	})

	t.Run("nil-input", func(t *testing.T) {
		out := NewCiphertext(params, 1, 5)
		outBefore := out.CopyNew()
		require.Error(t, eval.Rotate(nil, 1, out))
		require.True(t, out.Equal(outBefore))
		_, err := eval.RotateNew(nil, 1)
		require.Error(t, err)
	})
}

func TestPublicFastZeroSecretRotateRejectsConjugateInvariantRing(t *testing.T) {
	params := fastRotateTestParameters(t, ring.ConjugateInvariant)
	eval := NewEvaluator(params, nil)
	input, out := NewCiphertext(params, 1, 1), NewCiphertext(params, 1, 1)
	inputBefore, outBefore := input.CopyNew(), out.CopyNew()
	require.Error(t, eval.Rotate(input, 1, out))
	require.True(t, input.Equal(inputBefore))
	require.True(t, out.Equal(outBefore))
}
