package fast

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type trackingEvaluationKeySet struct {
	rlwe.EvaluationKeySet
	relinearizationKeyLookups int
}

func (keySet *trackingEvaluationKeySet) GetRelinearizationKey() (*rlwe.RelinearizationKey, error) {
	keySet.relinearizationKeyLookups++
	return keySet.EvaluationKeySet.GetRelinearizationKey()
}

func newPublicZeroSecretCiphertext(t *testing.T, params ckks.Parameters, encoder *ckks.Encoder, encryptor *rlwe.Encryptor, level int, values []complex128, scale float64) *rlwe.Ciphertext {
	t.Helper()
	pt := ckks.NewPlaintext(params, level)
	pt.Scale = rlwe.NewScale(scale)
	pt.IsNTT = true
	pt.LogDimensions = ring.Dimensions{Cols: params.LogMaxSlots()}
	require.NoError(t, encoder.Encode(values, pt))

	ct, err := encryptor.EncryptNew(pt)
	require.NoError(t, err)
	require.Equal(t, 1, ct.Degree())
	require.Equal(t, level, ct.Level())
	require.True(t, ct.IsNTT)
	require.False(t, ct.IsMontgomery)
	for q := 0; q <= level; q++ {
		for _, coefficient := range ct.Value[1].Coeffs[q] {
			require.Zero(t, coefficient, "Fast EncryptNew must produce zero c1 at active q%d", q)
		}
	}
	return ct
}

func requireDecodedProduct(t *testing.T, params ckks.Parameters, encoder *ckks.Encoder, decryptor *rlwe.Decryptor, got *rlwe.Ciphertext, expected []complex128, level int, expectedScale rlwe.Scale) {
	t.Helper()
	require.Equal(t, 1, got.Degree())
	require.Equal(t, level, got.Level())
	require.True(t, got.Scale.Equal(expectedScale))
	require.True(t, got.IsNTT)
	require.False(t, got.IsMontgomery)
	require.Equal(t, ring.Dimensions{Cols: params.LogMaxSlots()}, got.LogDimensions)
	require.Len(t, got.Value[0].Coeffs, level+1, "c0 must retain every active logical Q row")
	require.Len(t, got.Value[1].Coeffs, level+1, "c1 must retain every active logical Q row")
	for q := 0; q <= level; q++ {
		require.Len(t, got.Value[0].Coeffs[q], params.N())
		require.Len(t, got.Value[1].Coeffs[q], params.N())
	}
	for q := 0; q <= got.Level(); q++ {
		for _, coefficient := range got.Value[1].Coeffs[q] {
			require.Zero(t, coefficient, "keyless result c1 must be zero at active q%d", q)
		}
	}
	decoded := make([]complex128, len(expected))
	require.NoError(t, encoder.Decode(decryptor.DecryptNew(got), decoded))
	for i := range expected {
		require.LessOrEqual(t, cmplx.Abs(decoded[i]-expected[i]), 1e-4, "slot %d", i)
	}
}

func TestPublicCKKSEvaluatorKeylessZeroSecretMulRelin(t *testing.T) {
	params := testFastCKKSParameters(t)
	encoder := ckks.NewEncoder(params)
	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	decryptor := rlwe.NewDecryptor(params, sk)
	// A nil evaluation-key set is intentional: success proves MulRelin did not
	// access a RelinearizationKey or enter native GadgetProduct.
	eval := ckks.NewEvaluator(params, nil)
	// The supported keyless branch remains available when a normal CKKS
	// RelinearizationKey is present as well.
	standardRelinKey := kgen.GenRelinearizationKeyNew(sk)
	require.Equal(t, rlwe.KeyLayoutStandard, standardRelinKey.Layout)
	evalWithStandardKey := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(standardRelinKey))
	valuesA := []complex128{0.125 + 0.0625i, -0.25 + 0.125i, 0.0625 - 0.1875i, 0.3125 + 0i, -0.125 - 0.0625i, 0.1875 + 0.25i, -0.0625 + 0i, 0.25 - 0.125i}
	valuesB := []complex128{0.25 - 0.0625i, 0.125 + 0.1875i, -0.25 + 0.0625i, 0.0625 + 0.125i, 0.1875 + 0i, -0.125 + 0.0625i, 0.3125 - 0.1875i, -0.0625 - 0.125i}
	product := make([]complex128, len(valuesA))
	for i := range product {
		product[i] = valuesA[i] * valuesB[i]
	}

	for _, level := range []int{0, 3} {
		t.Run(map[int]string{0: "level0", 3: "level3"}[level], func(t *testing.T) {
			a := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, level, valuesA, math.Exp2(16))
			b := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, level, valuesB, math.Exp2(16))
			aBefore, bBefore := a.CopyNew(), b.CopyNew()

			got, err := eval.MulRelinNew(a, b)
			require.NoError(t, err)
			requireDecodedProduct(t, params, encoder, decryptor, got, product, level, a.Scale.Mul(b.Scale))
			gotWithKey, err := evalWithStandardKey.MulRelinNew(a, b)
			require.NoError(t, err)
			requireDecodedProduct(t, params, encoder, decryptor, gotWithKey, product, level, a.Scale.Mul(b.Scale))
			require.True(t, a.Equal(aBefore), "distinct-input MulRelinNew modified op0")
			require.True(t, b.Equal(bBefore), "distinct-input MulRelinNew modified op1")

			out := ckks.NewCiphertext(params, 1, level)
			require.NoError(t, eval.MulRelin(a, b, out))
			requireDecodedProduct(t, params, encoder, decryptor, out, product, level, a.Scale.Mul(b.Scale))

			square, err := eval.MulRelinNew(a, a)
			require.NoError(t, err)
			squareOracle := make([]complex128, len(valuesA))
			for i := range squareOracle {
				squareOracle[i] = valuesA[i] * valuesA[i]
			}
			requireDecodedProduct(t, params, encoder, decryptor, square, squareOracle, level, a.Scale.Mul(a.Scale))

			aAlias := a.CopyNew()
			require.NoError(t, eval.MulRelin(aAlias, b, aAlias))
			requireDecodedProduct(t, params, encoder, decryptor, aAlias, product, level, a.Scale.Mul(b.Scale))
			bAlias := b.CopyNew()
			require.NoError(t, eval.MulRelin(a, bAlias, bAlias))
			requireDecodedProduct(t, params, encoder, decryptor, bAlias, product, level, a.Scale.Mul(b.Scale))
		})
	}
}

func TestPublicCKKSEvaluatorKeylessZeroSecretMulRelinMetadata(t *testing.T) {
	params := testFastCKKSParameters(t)
	encoder := ckks.NewEncoder(params)
	sk := ckks.NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	eval := ckks.NewEvaluator(params, nil)
	values := []complex128{0.125, 0.25, -0.125, 0.0625, 0.1875, -0.0625, 0.25, -0.1875}
	a := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, 3, values, math.Exp2(16))
	b := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, 2, values, math.Exp2(15))
	a.LogDimensions = ring.Dimensions{Rows: 0, Cols: 2}
	b.LogDimensions = ring.Dimensions{Rows: 0, Cols: 3}
	out := ckks.NewCiphertext(params, 1, 3)

	require.NoError(t, eval.MulRelin(a, b, out))
	require.Equal(t, 2, out.Level())
	require.True(t, out.Scale.Equal(a.Scale.Mul(b.Scale)))
	require.Equal(t, ring.Dimensions{Rows: 0, Cols: 3}, out.LogDimensions)
	require.True(t, out.IsNTT)
	require.False(t, out.IsMontgomery)
}

func requireKeylessMulRelinRejectedWithoutMutation(t *testing.T, eval *ckks.Evaluator, a, b, out *rlwe.Ciphertext) {
	t.Helper()
	aBefore, bBefore, outBefore := a.CopyNew(), b.CopyNew(), out.CopyNew()
	require.Error(t, eval.MulRelin(a, b, out))
	requireCiphertextByteIdentical(t, a, aBefore, "rejected MulRelin modified op0")
	requireCiphertextByteIdentical(t, b, bBefore, "rejected MulRelin modified op1")
	requireCiphertextByteIdentical(t, out, outBefore, "rejected MulRelin modified output")
}

func requireCiphertextByteIdentical(t *testing.T, got, want *rlwe.Ciphertext, message string) {
	t.Helper()
	gotBytes, err := got.MarshalBinary()
	require.NoError(t, err)
	wantBytes, err := want.MarshalBinary()
	require.NoError(t, err)
	require.Equal(t, wantBytes, gotBytes, message)
}

func TestPublicCKKSEvaluatorKeylessZeroSecretMulRelinRejectsUnsafeInputs(t *testing.T) {
	params := testFastCKKSParameters(t)
	encoder := ckks.NewEncoder(params)
	sk := ckks.NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	eval := ckks.NewEvaluator(params, nil)
	values := []complex128{0.125, 0.25, -0.125, 0.0625, 0.1875, -0.0625, 0.25, -0.1875}
	newPair := func() (*rlwe.Ciphertext, *rlwe.Ciphertext, *rlwe.Ciphertext) {
		a := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, 3, values, math.Exp2(16))
		b := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, 3, values, math.Exp2(16))
		out := ckks.NewCiphertext(params, 1, 3)
		out.Scale = rlwe.NewScale(17)
		return a, b, out
	}

	t.Run("nonzero c1", func(t *testing.T) {
		a, b, out := newPair()
		a.Value[1].Coeffs[1][0] = 1
		aBefore, bBefore, outBefore := a.CopyNew(), b.CopyNew(), out.CopyNew()
		err := eval.MulRelin(a, b, out)
		require.ErrorContains(t, err, "nonzero c1")
		requireCiphertextByteIdentical(t, a, aBefore, "nonzero-c1 rejection modified op0")
		requireCiphertextByteIdentical(t, b, bBefore, "nonzero-c1 rejection modified op1")
		requireCiphertextByteIdentical(t, out, outBefore, "nonzero-c1 rejection modified output")
	})

	t.Run("compact active Q backing", func(t *testing.T) {
		a, b, out := newPair()
		a.Value[0].Coeffs[3] = nil
		requireKeylessMulRelinRejectedWithoutMutation(t, eval, a, b, out)
	})

	t.Run("logical level beyond configured Q chain", func(t *testing.T) {
		a, b, out := newPair()
		for component := range a.Value {
			a.Value[component].Coeffs = append(a.Value[component].Coeffs, make([]uint64, params.N()))
		}
		require.Equal(t, 4, a.Level())
		requireKeylessMulRelinRejectedWithoutMutation(t, eval, a, b, out)
	})

	t.Run("unsupported Montgomery representation", func(t *testing.T) {
		a, b, out := newPair()
		params.RingQ().AtLevel(3).MForm(a.Value[0], a.Value[0])
		params.RingQ().AtLevel(3).MForm(b.Value[0], b.Value[0])
		a.IsMontgomery, b.IsMontgomery = true, true
		requireKeylessMulRelinRejectedWithoutMutation(t, eval, a, b, out)
	})

	t.Run("unsupported non-NTT representation", func(t *testing.T) {
		a, b, out := newPair()
		a.IsNTT, b.IsNTT = false, false
		requireKeylessMulRelinRejectedWithoutMutation(t, eval, a, b, out)
	})

	t.Run("degree greater than one", func(t *testing.T) {
		degreeTwo := ckks.NewCiphertext(params, 2, 3)
		degreeOne := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, 3, values, math.Exp2(16))
		_, err := eval.MulRelinNew(degreeTwo, degreeOne)
		require.Error(t, err)
	})
}

func TestPublicCKKSEvaluatorFastMulRelinFailsClosedWithStandardRelinKey(t *testing.T) {
	params := testFastCKKSParameters(t)
	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	relinKey := kgen.GenRelinearizationKeyNew(sk)
	require.Equal(t, rlwe.KeyLayoutStandard, relinKey.Layout)
	keySet := &trackingEvaluationKeySet{EvaluationKeySet: rlwe.NewMemEvaluationKeySet(relinKey)}
	eval := ckks.NewEvaluator(params, keySet)
	encoder := ckks.NewEncoder(params)
	encryptor := rlwe.NewEncryptor(params, sk)
	values := []complex128{0.125, 0.25, -0.125, 0.0625, 0.1875, -0.0625, 0.25, -0.1875}
	newPair := func() (*rlwe.Ciphertext, *rlwe.Ciphertext, *rlwe.Ciphertext) {
		a := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, 3, values, math.Exp2(16))
		b := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, 3, values, math.Exp2(16))
		out := ckks.NewCiphertext(params, 1, 3)
		out.Scale = rlwe.NewScale(17)
		out.LogDimensions = ring.Dimensions{Rows: 1, Cols: 2}
		return a, b, out
	}

	t.Run("nonzero c1 rejects despite available native key", func(t *testing.T) {
		a, b, out := newPair()
		a.Value[1].Coeffs[1][0] = 1
		aBefore, bBefore, outBefore := a.CopyNew(), b.CopyNew(), out.CopyNew()
		err := eval.MulRelin(a, b, out)
		require.ErrorContains(t, err, "op0 has nonzero c1 at q1")
		requireCiphertextByteIdentical(t, a, aBefore, "rejected MulRelin modified op0")
		requireCiphertextByteIdentical(t, b, bBefore, "rejected MulRelin modified op1")
		requireCiphertextByteIdentical(t, out, outBefore, "rejected MulRelin modified output metadata or data")

		_, err = eval.MulRelinNew(a, b)
		require.ErrorContains(t, err, "op0 has nonzero c1 at q1")
		requireCiphertextByteIdentical(t, a, aBefore, "rejected MulRelinNew modified op0")
		requireCiphertextByteIdentical(t, b, bBefore, "rejected MulRelinNew modified op1")
		require.Zero(t, keySet.relinearizationKeyLookups, "rejected Fast input must not inspect the available native relin key")
	})

	t.Run("higher degree rejects despite available native key", func(t *testing.T) {
		degreeTwo := ckks.NewCiphertext(params, 2, 3)
		degreeOne := newPublicZeroSecretCiphertext(t, params, encoder, encryptor, 3, values, math.Exp2(16))
		out := ckks.NewCiphertext(params, 1, 3)
		out.Scale = rlwe.NewScale(19)
		degreeTwoBefore, degreeOneBefore, outBefore := degreeTwo.CopyNew(), degreeOne.CopyNew(), out.CopyNew()
		err := eval.MulRelin(degreeTwo, degreeOne, out)
		require.ErrorContains(t, err, "ciphertext inputs must both have degree 1")
		requireCiphertextByteIdentical(t, degreeTwo, degreeTwoBefore, "degree rejection modified op0")
		requireCiphertextByteIdentical(t, degreeOne, degreeOneBefore, "degree rejection modified op1")
		requireCiphertextByteIdentical(t, out, outBefore, "degree rejection modified output")

		_, err = eval.MulRelinNew(degreeTwo, degreeOne)
		require.ErrorContains(t, err, "ciphertext inputs must both have degree 1")
		requireCiphertextByteIdentical(t, degreeTwo, degreeTwoBefore, "degree rejection through MulRelinNew modified op0")
		requireCiphertextByteIdentical(t, degreeOne, degreeOneBefore, "degree rejection through MulRelinNew modified op1")
		require.Zero(t, keySet.relinearizationKeyLookups, "rejected degree must not inspect the available native relin key")
	})

	t.Run("nil ciphertexts return errors without panic", func(t *testing.T) {
		a, b, out := newPair()
		_, err := eval.MulRelinNew(nil, b)
		require.ErrorContains(t, err, "ciphertexts cannot be nil")

		var nilCiphertext *rlwe.Ciphertext
		_, err = eval.MulRelinNew(a, nilCiphertext)
		require.ErrorContains(t, err, "ciphertexts cannot be nil")

		err = eval.MulRelin(nil, b, out)
		require.ErrorContains(t, err, "ciphertexts cannot be nil")
	})
}
