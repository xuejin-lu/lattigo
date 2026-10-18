package rlwe

import (
	"math"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

type aZeroSamplerCalls struct {
	read       int
	readAndAdd int
	samples    []ring.Poly
}

type aZeroCountingSampler struct {
	ring.Sampler
	calls *aZeroSamplerCalls
}

func (s aZeroCountingSampler) Read(pol ring.Poly) {
	s.calls.read++
	s.Sampler.Read(pol)
	s.calls.samples = append(s.calls.samples, *pol.CopyNew())
}

func (s aZeroCountingSampler) ReadNew() ring.Poly {
	s.calls.read++
	return s.Sampler.ReadNew()
}

func (s aZeroCountingSampler) ReadAndAdd(pol ring.Poly) {
	s.calls.readAndAdd++
	s.Sampler.ReadAndAdd(pol)
}

func (s aZeroCountingSampler) AtLevel(level int) ring.Sampler {
	return aZeroCountingSampler{Sampler: s.Sampler.AtLevel(level), calls: s.calls}
}

func aZeroParameters(t *testing.T, withP bool) Parameters {
	t.Helper()
	literal := testInsecure[0].ParametersLiteral
	if !withP {
		literal.P = nil
	}
	params, err := NewParametersFromLiteral(literal)
	require.NoError(t, err)
	return params
}

func aZeroKeyPair(t *testing.T, params Parameters, seed []byte) (*KeyGenerator, *SecretKey, *PublicKey) {
	t.Helper()
	prng, err := sampling.NewKeyedPRNG(seed)
	require.NoError(t, err)
	kgen := NewKeyGenerator(params)
	kgen.Encryptor = newTestEncryptorWithKeyedPRNG(params, nil, prng)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	return kgen, sk, pk
}

func requireAZeroPoly(t *testing.T, pol ring.Poly) {
	t.Helper()
	for level := range pol.Coeffs {
		for _, coefficient := range pol.Coeffs[level] {
			require.Zero(t, coefficient)
		}
	}
}

func requireANonZeroPoly(t *testing.T, pol ring.Poly) {
	t.Helper()
	for level := range pol.Coeffs {
		for _, coefficient := range pol.Coeffs[level] {
			if coefficient != 0 {
				return
			}
		}
	}
	t.Fatal("expected polynomial to contain non-zero coefficients")
}

func requireAZeroQP(t *testing.T, pol ringqp.Poly) {
	t.Helper()
	requireAZeroPoly(t, pol.Q)
	requireAZeroPoly(t, pol.P)
}

func requireANonZeroQP(t *testing.T, pol ringqp.Poly) {
	t.Helper()
	if hasNonZeroCoefficient(pol.Q) || hasNonZeroCoefficient(pol.P) {
		return
	}
	t.Fatal("expected QP polynomial to contain non-zero coefficients")
}

func hasNonZeroCoefficient(pol ring.Poly) bool {
	for level := range pol.Coeffs {
		for _, coefficient := range pol.Coeffs[level] {
			if coefficient != 0 {
				return true
			}
		}
	}
	return false
}

func TestAZeroPreservesSecretAndEncryptionNoise(t *testing.T) {
	params := aZeroParameters(t, true)
	seed := []byte("standard-a-zero-key-seed")

	// The secret-key output still exactly matches the unchanged Xs sampler.
	actualPRNG, err := sampling.NewKeyedPRNG(seed)
	require.NoError(t, err)
	actualKgen := NewKeyGenerator(params)
	actualKgen.Encryptor = newTestEncryptorWithKeyedPRNG(params, nil, actualPRNG)
	sk := actualKgen.GenSecretKeyNew()

	referencePRNG, err := sampling.NewKeyedPRNG(seed)
	require.NoError(t, err)
	referenceEnc := newTestEncryptorWithKeyedPRNG(params, nil, referencePRNG)
	expectedSk := NewSecretKey(params)
	referenceEnc.xsSampler.AtLevel(params.MaxLevelQ()).Read(expectedSk.Value.Q)
	if params.MaxLevelP() >= 0 {
		params.RingQP().AtLevel(params.MaxLevelQ(), params.MaxLevelP()).ExtendBasisSmallNormAndCenter(
			expectedSk.Value.Q, params.MaxLevelP(), expectedSk.Value.Q, expectedSk.Value.P)
	}
	params.RingQP().NTT(expectedSk.Value, expectedSk.Value)
	params.RingQP().MForm(expectedSk.Value, expectedSk.Value)
	require.True(t, sk.Value.Equal(&expectedSk.Value), "secret-key Xs sampling changed")
	requireANonZeroQP(t, sk.Value)

	// Replay the original a draw before e. Equality with this deterministic e
	// proves both the draw order and the public-key noise path were preserved.
	_, publicKeySk, pk := aZeroKeyPair(t, params, seed)
	requireAZeroPoly(t, pk.Value[1].Q)
	requireAZeroPoly(t, pk.Value[1].P)

	replayPRNG, err := sampling.NewKeyedPRNG(seed)
	require.NoError(t, err)
	replayEnc := newTestEncryptorWithKeyedPRNG(params, nil, replayPRNG)
	replayKgen := NewKeyGenerator(params)
	replayKgen.Encryptor = replayEnc
	replayKgen.GenSecretKeyNew()
	sampledA := params.RingQP().AtLevel(params.MaxLevelQ(), params.MaxLevelP()).NewPoly()
	replayEnc.uniformSampler.AtLevel(params.MaxLevelQ(), params.MaxLevelP()).Read(sampledA)
	requireANonZeroQP(t, sampledA)
	expectedPkError := params.RingQP().AtLevel(params.MaxLevelQ(), params.MaxLevelP()).NewPoly()
	replayEnc.xeSampler.AtLevel(params.MaxLevelQ()).Read(expectedPkError.Q)
	params.RingQP().AtLevel(params.MaxLevelQ(), params.MaxLevelP()).ExtendBasisSmallNormAndCenter(
		expectedPkError.Q, params.MaxLevelP(), expectedPkError.Q, expectedPkError.P)
	params.RingQP().AtLevel(params.MaxLevelQ(), params.MaxLevelP()).NTT(expectedPkError, expectedPkError)
	params.RingQP().AtLevel(params.MaxLevelQ(), params.MaxLevelP()).MForm(expectedPkError, expectedPkError)
	require.True(t, pk.Value[0].Equal(&expectedPkError), "public-key error draw or PRNG order changed")
	require.False(t, math.IsNaN(NoisePublicKey(pk, publicKeySk, params)))
	require.False(t, math.IsInf(NoisePublicKey(pk, publicKeySk, params), 0))
	requireANonZeroQP(t, pk.Value[0])

	// Secret-key encryption consumes its original Q-only uniform draw, then
	// retains the following xe sample as c0 while clearing only c1=a.
	secretEncrypt := func(seed []byte) (*Ciphertext, *aZeroSamplerCalls) {
		prng, err := sampling.NewKeyedPRNG(seed)
		require.NoError(t, err)
		enc := newTestEncryptorWithKeyedPRNG(params, sk, prng)
		calls := &aZeroSamplerCalls{}
		enc.xeSampler = aZeroCountingSampler{Sampler: enc.xeSampler, calls: calls}
		ct := NewCiphertext(params, 1, params.MaxLevelQ())
		require.NoError(t, enc.EncryptZero(ct))
		return ct, calls
	}

	encSeed := []byte("standard-a-zero-secret-encryption")
	ct, calls := secretEncrypt(encSeed)
	require.Equal(t, 1, calls.read)
	require.Zero(t, calls.readAndAdd)
	requireAZeroPoly(t, ct.Value[1])
	requireANonZeroPoly(t, ct.Value[0])

	replayEncPRNG, err := sampling.NewKeyedPRNG(encSeed)
	require.NoError(t, err)
	replaySecretEnc := newTestEncryptorWithKeyedPRNG(params, sk, replayEncPRNG)
	sampledAQ := params.RingQ().AtLevel(params.MaxLevelQ()).NewPoly()
	replaySecretEnc.uniformSampler.AtLevel(params.MaxLevelQ(), -1).Read(ringqp.Poly{Q: sampledAQ})
	requireANonZeroPoly(t, sampledAQ)
	expectedSecretError := params.RingQ().AtLevel(params.MaxLevelQ()).NewPoly()
	replaySecretEnc.xeSampler.AtLevel(params.MaxLevelQ()).Read(expectedSecretError)
	params.RingQ().AtLevel(params.MaxLevelQ()).NTT(expectedSecretError, expectedSecretError)
	require.True(t, params.RingQ().Equal(ct.Value[0], expectedSecretError), "secret-key error draw or PRNG order changed")

	repeatedCT, _ := secretEncrypt(encSeed)
	require.True(t, ct.Equal(repeatedCT), "fixed-seed secret-key encryption was not reproducible")
	noiseStd, _, _ := Norm(ct, NewDecryptor(params, sk))
	require.False(t, math.IsNaN(noiseStd))
	require.False(t, math.IsInf(noiseStd, 0))
}

func TestAZeroPreservesPublicEncryptionErrorDraws(t *testing.T) {
	for _, withP := range []bool{true, false} {
		name := "Q-only"
		if withP {
			name = "QP"
		}
		t.Run(name, func(t *testing.T) {
			params := aZeroParameters(t, withP)
			_, _, pk := aZeroKeyPair(t, params, []byte("standard-a-zero-pk-"+name))
			requireAZeroPoly(t, pk.Value[1].Q)
			requireAZeroPoly(t, pk.Value[1].P)

			encrypt := func() (*Ciphertext, *aZeroSamplerCalls) {
				prng, err := sampling.NewKeyedPRNG([]byte("standard-a-zero-public-encryption-" + name))
				require.NoError(t, err)
				enc := newTestEncryptorWithKeyedPRNG(params, pk, prng)
				calls := &aZeroSamplerCalls{}
				enc.xeSampler = aZeroCountingSampler{Sampler: enc.xeSampler, calls: calls}
				ct := NewCiphertext(params, 1, params.MaxLevelQ())
				require.NoError(t, enc.EncryptZero(ct))
				return ct, calls
			}

			ct, calls := encrypt()
			require.Equal(t, 2, calls.read+calls.readAndAdd, "public-key e0/e1 must both be sampled")
			require.Len(t, calls.samples, 2)
			requireANonZeroPoly(t, calls.samples[0])
			requireANonZeroPoly(t, calls.samples[1])
			if !withP {
				requireANonZeroPoly(t, ct.Value[0])
				requireANonZeroPoly(t, ct.Value[1]) // In Q-only mode c1 is the preserved e1.
			}

			repeatedCT, repeatedCalls := encrypt()
			require.Equal(t, 2, repeatedCalls.read+repeatedCalls.readAndAdd)
			require.True(t, ct.Equal(repeatedCT), "fixed-seed public-key encryption was not reproducible")
			require.Len(t, repeatedCalls.samples, 2)
			for i := range calls.samples {
				require.True(t, calls.samples[i].Equal(&repeatedCalls.samples[i]), "fixed-seed e%d draw was not reproducible", i)
			}
		})
	}
}

func requireAZeroEvaluationKey(t *testing.T, evk *EvaluationKey) {
	t.Helper()
	require.Equal(t, 1, evk.Degree())
	for i := range evk.Value {
		for j := range evk.Value[i] {
			requireANonZeroQP(t, evk.Value[i][j][0])
			requireAZeroQP(t, evk.Value[i][j][1])
		}
	}
}

func requireAZeroEvaluationKeyNoise(t *testing.T, noise float64, params Parameters, evk *EvaluationKey) {
	t.Helper()
	require.False(t, math.IsNaN(noise))
	require.False(t, math.IsInf(noise, 0))
	noiseBound := math.Log2(math.Sqrt(float64(len(evk.Value)))*params.NoiseFreshSK()) + 2
	require.LessOrEqual(t, noise, noiseBound, "evaluation-key first component lost its gadget payload or has excessive noise")
}

func TestAZeroEvaluationKeysAndCompressedExpansion(t *testing.T) {
	params := aZeroParameters(t, true)
	keySeed := []byte("standard-a-zero-evaluation-keys")
	kgen, sk, _ := aZeroKeyPair(t, params, keySeed)
	skOut := kgen.GenSecretKeyNew()

	evk := kgen.GenEvaluationKeyNew(sk, skOut)
	requireAZeroEvaluationKey(t, evk)
	evalKeyNoise := NoiseEvaluationKey(evk, sk, skOut, params)
	requireAZeroEvaluationKeyNoise(t, evalKeyNoise, params, evk)

	repeatedKgen, repeatedSk, _ := aZeroKeyPair(t, params, keySeed)
	repeatedSkOut := repeatedKgen.GenSecretKeyNew()
	repeatedEvk := repeatedKgen.GenEvaluationKeyNew(repeatedSk, repeatedSkOut)
	repeatedEvkBytes, err := repeatedEvk.MarshalBinary()
	require.NoError(t, err)
	evkBytes, err := evk.MarshalBinary()
	require.NoError(t, err)
	require.Equal(t, evkBytes, repeatedEvkBytes, "fixed-seed evaluation-key error sampling was not reproducible")
	require.Equal(t, evalKeyNoise, NoiseEvaluationKey(repeatedEvk, repeatedSk, repeatedSkOut, params))

	rlk := kgen.GenRelinearizationKeyNew(sk)
	requireAZeroEvaluationKey(t, &rlk.EvaluationKey)
	requireAZeroEvaluationKeyNoise(t, NoiseRelinearizationKey(rlk, sk, params), params, &rlk.EvaluationKey)

	gk := kgen.GenGaloisKeyNew(params.GaloisElement(3), sk)
	requireAZeroEvaluationKey(t, &gk.EvaluationKey)
	requireAZeroEvaluationKeyNoise(t, NoiseGaloisKey(gk, sk, params), params, &gk.EvaluationKey)

	// Dimension-switching and standard/conjugate-invariant ring-swap keys use
	// the same generic EncryptZero-based gadget path.
	paramsSmall, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: params.LogN() - 1, Q: params.Q(), P: params.P(), NTTFlag: params.NTTFlag(),
	})
	require.NoError(t, err)
	skSmall := NewKeyGenerator(paramsSmall).GenSecretKeyNew()
	for _, dimensionKey := range []*EvaluationKey{
		kgen.GenEvaluationKeyNew(skSmall, sk),
		kgen.GenEvaluationKeyNew(sk, skSmall),
	} {
		requireAZeroEvaluationKey(t, dimensionKey)
	}

	paramsCI, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: params.LogN() - 1, Q: params.Q(), P: params.P(), NTTFlag: params.NTTFlag(), RingType: ring.ConjugateInvariant,
	})
	require.NoError(t, err)
	skCI := NewKeyGenerator(paramsCI).GenSecretKeyNew()
	stdToCI, ciToStd := kgen.GenEvaluationKeysForRingSwapNew(sk, skCI)
	requireAZeroEvaluationKey(t, stdToCI)
	requireAZeroEvaluationKey(t, ciToStd)

	compressed := kgen.GenEvaluationKeyNew(sk, skOut, EvaluationKeyParameters{Compressed: true})
	require.True(t, compressed.IsCompressed())
	firstBeforeExpand := make([][]ringqp.Poly, len(compressed.Value))
	for i := range compressed.Value {
		firstBeforeExpand[i] = make([]ringqp.Poly, len(compressed.Value[i]))
		for j := range compressed.Value[i] {
			firstBeforeExpand[i][j] = *compressed.Value[i][j][0].CopyNew()
			requireANonZeroQP(t, firstBeforeExpand[i][j])
		}
	}
	require.NoError(t, compressed.Expand(params, nil))
	require.False(t, compressed.IsCompressed())
	requireAZeroEvaluationKey(t, compressed)
	for i := range compressed.Value {
		for j := range compressed.Value[i] {
			require.True(t, firstBeforeExpand[i][j].Equal(&compressed.Value[i][j][0]), "compressed expansion changed the first component")
		}
	}
	requireAZeroEvaluationKeyNoise(t, NoiseEvaluationKey(compressed, sk, skOut, params), params, compressed)

	// A non-zero first component and finite, bounded residual noise demonstrate
	// that key payload and e survive while its Q/P second component remains 0.
	for keyIndex, key := range []*EvaluationKey{evk, &rlk.EvaluationKey, &gk.EvaluationKey, compressed} {
		for i := range key.Value {
			for j := range key.Value[i] {
				require.True(t, hasNonZeroCoefficient(key.Value[i][j][0].Q) || hasNonZeroCoefficient(key.Value[i][j][0].P),
					"key %d entry (%d,%d) lost its first component", keyIndex, i, j)
			}
		}
	}

	// Public-key encryption keeps c1=e1 non-zero, allowing the existing
	// evaluation-key application path to exercise a real key switch.
	paramsNoP := aZeroParameters(t, false)
	kgenNoP, skNoP, pkNoP := aZeroKeyPair(t, paramsNoP, []byte("standard-a-zero-keyswitch-pk"))
	skOutNoP := kgenNoP.GenSecretKeyNew()
	evkNoP := kgenNoP.GenEvaluationKeyNew(skNoP, skOutNoP)
	enc := NewEncryptor(paramsNoP, pkNoP)
	ct := NewCiphertext(paramsNoP, 1, paramsNoP.MaxLevelQ())
	require.NoError(t, enc.EncryptZero(ct))
	requireANonZeroPoly(t, ct.Value[1])
	eval := NewEvaluator(paramsNoP, nil)
	require.NoError(t, eval.ApplyEvaluationKey(ct, evkNoP, ct))
	noiseStd, _, _ := Norm(ct, NewDecryptor(paramsNoP, skOutNoP))
	require.False(t, math.IsNaN(noiseStd))
	require.False(t, math.IsInf(noiseStd, 0))
}
