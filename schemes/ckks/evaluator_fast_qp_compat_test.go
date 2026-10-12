package ckks

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

func TestGenericFastCKKSKeyGeneratorProducesZeroSecretEvaluationKeys(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	original := sk.CopyNew()

	rlk := kgen.GenRelinearizationKeyNew(sk)
	galEl := params.GaloisElementForRotation(1)
	gk := kgen.GenGaloisKeyNew(galEl, sk)
	evaluationKey := kgen.GenEvaluationKeyNew(sk, sk)

	for name, key := range map[string]*rlwe.EvaluationKey{
		"relinearization": &rlk.EvaluationKey,
		"galois":          &gk.EvaluationKey,
		"evaluation":      evaluationKey,
	} {
		require.Equal(t, rlwe.KeyLayoutFast, key.Layout, name)
		require.Equal(t, rlwe.KeyLayoutFast, key.GadgetCiphertext.Layout, name)
		for i := range key.Value {
			for j := range key.Value[i] {
				require.True(t, fastQPolyIsZero(key.Value[i][j][1].Q, len(key.Value[i][j][1].Q.Coeffs)), name)
				require.True(t, fastQPolyIsZero(key.Value[i][j][1].P, len(key.Value[i][j][1].P.Coeffs)), name)
			}
		}
	}
	require.Equal(t, original, sk, "generic Fast key generation must not mutate the caller's secret key")
}

func TestFastQPCompatibilityDecomposeAutomorphismAndModDown(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	galEl := params.GaloisElementForRotation(1)
	gk := kgen.GenGaloisKeyNew(galEl, sk)
	eval := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
	encoder := NewEncoder(params)
	encryptor := rlwe.NewEncryptor(params, sk)
	level := 3
	input := fastRotateTestEncrypt(t, params, encoder, encryptor, fastRotateTestValues(4), 4, level, true)

	beta := params.BaseRNSDecompositionVectorSize(level, gk.LevelP())
	decomposition := make([]ringqp.Poly, beta)
	for i := range decomposition {
		decomposition[i] = params.RingQP().AtLevel(level, gk.LevelP()).NewPoly()
	}
	eval.DecomposeNTT(level, gk.LevelP(), gk.LevelP()+1, input.Value[1], input.IsNTT, decomposition)
	require.True(t, fastQPDecompositionIsZero(decomposition, level, gk.LevelP()))
	require.Equal(t, uint64(1), eval.fastQPCompat.counterSnapshot().decomposeElisions)

	qp := rlwe.NewElementExtended(params, 1, level, gk.LevelP())
	require.NoError(t, eval.AutomorphismHoistedLazy(level, input, decomposition, galEl, qp))
	counters := eval.fastQPCompat.counterSnapshot()
	require.Equal(t, uint64(1), counters.automorphisms, "the compatibility path must perform a Galois automorphism")
	require.Equal(t, uint64(1), counters.keySwitchElisions, "zero c1 must elide key switching")
	ringQ := params.RingQ().AtLevel(level)
	rotated := ringQ.NewPoly()
	scaled := ringQ.NewPoly()
	ringQ.MulScalarBigint(input.Value[0], params.RingP().ModulusAtLevel[gk.LevelP()], scaled)
	ringQ.AutomorphismNTTWithIndex(scaled, eval.AutomorphismIndex(galEl), rotated)
	for row := 0; row <= level; row++ {
		require.Equal(t, rotated.Coeffs[row], qp.Value[0].Q.Coeffs[row], "P-scaled Q output row %d", row)
		require.Zero(t, qp.Value[0].P.Coeffs[0][0])
		for _, value := range qp.Value[1].Q.Coeffs[row] {
			require.Zero(t, value)
		}
	}
	for row := 0; row <= gk.LevelP(); row++ {
		for _, value := range qp.Value[0].P.Coeffs[row] {
			require.Zero(t, value)
		}
		for _, value := range qp.Value[1].P.Coeffs[row] {
			require.Zero(t, value)
		}
	}

	got := rlwe.NewCiphertext(params, 1, level)
	eval.ModDownQPtoQNTT(level, gk.LevelP(), qp.Value[0].Q, qp.Value[0].P, got.Value[0])
	eval.ModDownQPtoQNTT(level, gk.LevelP(), qp.Value[1].Q, qp.Value[1].P, got.Value[1])
	want := ringQ.NewPoly()
	ringQ.AutomorphismNTTWithIndex(input.Value[0], eval.AutomorphismIndex(galEl), want)
	for row := 0; row <= level; row++ {
		require.Equal(t, want.Coeffs[row], got.Value[0].Coeffs[row], "ModDown row %d", row)
	}
	require.True(t, fastQPolyIsZero(got.Value[1], level+1))
}

func TestPublicFastAutomorphismUsesPrefixCoreAndVerifiedFullQPath(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

	t.Run("ordinary-Conjugate-uses-shared-prefix-core", func(t *testing.T) {
		eval := NewEvaluator(params, nil)
		recorder := &recordingPublicAutomorphismCore{AutomorphismCore: eval.automorphismCore}
		eval.automorphismCore = recorder
		input := fastRotateTestEncrypt(t, params, encoder, rlwe.NewEncryptor(params, sk), fastRotateTestValues(4), 4, 3, true)
		output := rlwe.NewCiphertext(params, 1, input.Level())
		output.IsNTT = input.IsNTT

		require.NoError(t, eval.Conjugate(input, output))
		require.Equal(t, []int{4}, recorder.rows, "ordinary Conjugate must dispatch through the maintained Q-prefix core")
		require.Equal(t, []int{2}, recorder.pairCount)
		ringQ := params.RingQ().AtLevel(input.Level())
		want := ringQ.NewPoly()
		ringQ.AutomorphismNTT(input.Value[0], params.GaloisElementForComplexConjugation(), want)
		for row := 0; row <= input.Level(); row++ {
			require.Equal(t, want.Coeffs[row], output.Value[0].Coeffs[row])
			for _, coefficient := range output.Value[1].Coeffs[row] {
				require.Zero(t, coefficient)
			}
		}
		require.Equal(t, *input.MetaData, *output.MetaData)
	})

	t.Run("verified-full-Q-preserves-rows-and-rejects-incoherent-dormant-row", func(t *testing.T) {
		const level = 5
		eval := NewEvaluator(params, nil)
		input := rlwe.NewCiphertext(params, 1, level)
		input.IsNTT = false
		input.Scale = rlwe.NewScale(1 << fastRotateTestLogScale)
		// A small coefficient-domain polynomial has an unambiguous centered lift
		// into every Q row and is enough to exercise the actual Galois permutation.
		for row := 0; row <= level; row++ {
			input.Value[0].Coeffs[row][1] = 7
			input.Value[0].Coeffs[row][2] = 3
		}
		output := rlwe.NewCiphertext(params, 1, level)
		output.IsNTT = false
		galEl := params.GaloisElementForComplexConjugation()
		ringQ := params.RingQ().AtLevel(level)
		want := ringQ.NewPoly()
		ringQ.Automorphism(input.Value[0], galEl, want)

		require.NoError(t, eval.Conjugate(input, output))
		require.Equal(t, level, output.Level())
		require.False(t, output.IsNTT)
		require.True(t, fullFastPolyRows(output.Value[0], level, params.N()), "verified full-Q rows must remain materialized for BTS callers")
		for row := 0; row <= level; row++ {
			require.Equal(t, want.Coeffs[row], output.Value[0].Coeffs[row], "q%d", row)
			for _, coefficient := range output.Value[1].Coeffs[row] {
				require.Zero(t, coefficient)
			}
		}
		require.Equal(t, *input.MetaData, *output.MetaData)
		require.Equal(t, uint64(1), eval.fastQPCompat.counterSnapshot().automorphisms)

		bad := input.CopyNew()
		bad.Value[0].Coeffs[4][1]++
		badBefore := bad.CopyNew()
		outputBefore := output.CopyNew()
		require.ErrorContains(t, eval.Conjugate(bad, output), "complete Q input is not authorized")
		require.True(t, bad.Equal(badBefore))
		require.True(t, output.Equal(outputBefore), "failed full-Q verification must not modify output")
	})
}

func TestFastQPCompatibilityRejectsNonzeroInputsTransactionally(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	galEl := params.GaloisElementForRotation(1)
	gk := kgen.GenGaloisKeyNew(galEl, sk)
	eval := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
	encoder := NewEncoder(params)
	encryptor := rlwe.NewEncryptor(params, sk)
	input := fastRotateTestEncrypt(t, params, encoder, encryptor, fastRotateTestValues(4), 4, 3, true)

	inputBefore := input.CopyNew()
	input.Value[1].Coeffs[0][0] = 1
	mutatedBefore := input.CopyNew()
	beta := params.BaseRNSDecompositionVectorSize(input.Level(), gk.LevelP())
	decomposition := make([]ringqp.Poly, beta)
	for i := range decomposition {
		decomposition[i] = params.RingQP().AtLevel(input.Level(), gk.LevelP()).NewPoly()
		decomposition[i].Q.Coeffs[0][0] = 17
	}
	output := rlwe.NewElementExtended(params, 1, input.Level(), gk.LevelP())
	for component := range output.Value {
		for row := range output.Value[component].Q.Coeffs {
			for i := range output.Value[component].Q.Coeffs[row] {
				output.Value[component].Q.Coeffs[row][i] = 23
			}
		}
	}
	outputBefore := output.CopyNew()
	err := eval.AutomorphismHoistedLazy(input.Level(), input, decomposition, galEl, output)
	require.Error(t, err)
	require.Equal(t, outputBefore, output, "failed automorphism must not modify the existing QP output")
	require.Equal(t, mutatedBefore, input, "failed automorphism must not modify its input")
	require.NotEqual(t, inputBefore, input, "test setup must actually exercise nonzero c1")

	// DecomposeNTT cannot return an error by API contract. It must leave the
	// caller's buffers untouched and surface the rejection at the next API
	// boundary rather than panic or manufacture zero digits.
	for i := range decomposition {
		decomposition[i].Q.Zero()
		decomposition[i].Q.Coeffs[0][0] = 31
	}
	eval.DecomposeNTT(input.Level(), gk.LevelP(), gk.LevelP()+1, input.Value[1], input.IsNTT, decomposition)
	for i := range decomposition {
		require.Equal(t, uint64(31), decomposition[i].Q.Coeffs[0][0])
	}
	err = eval.Trace(input, params.LogMaxSlots(), input)
	require.Error(t, err)
	require.Equal(t, mutatedBefore, input)
	for i := range decomposition {
		decomposition[i].Q.Coeffs[0][0] = 47
	}
	eval.DecomposeNTT(input.Level(), gk.LevelP(), gk.LevelP()+1, ring.Poly{}, input.IsNTT, decomposition)
	for i := range decomposition {
		require.Equal(t, uint64(47), decomposition[i].Q.Coeffs[0][0], "malformed DecomposeNTT must leave output untouched")
	}
	require.Error(t, eval.Trace(input, params.LogMaxSlots(), input), "malformed void DecomposeNTT error must be surfaced")

}

func TestFastQPCompatibilityApplyEvaluationKeyAndGadgetProduct(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	evaluationKey := kgen.GenEvaluationKeyNew(sk, sk)
	eval := NewEvaluator(params, nil)
	encoder := NewEncoder(params)
	encryptor := rlwe.NewEncryptor(params, sk)
	input := fastRotateTestEncrypt(t, params, encoder, encryptor, fastRotateTestValues(4), 4, 3, true)

	output := rlwe.NewCiphertext(params, 1, input.Level())
	require.NoError(t, eval.ApplyEvaluationKey(input, evaluationKey, output))
	require.Equal(t, uint64(1), eval.fastQPCompat.counterSnapshot().applyEvaluationKeyElisions)
	for row := 0; row <= input.Level(); row++ {
		require.Equal(t, input.Value[0].Coeffs[row], output.Value[0].Coeffs[row])
		for _, value := range output.Value[1].Coeffs[row] {
			require.Zero(t, value)
		}
	}

	badInput := input.CopyNew()
	badInput.Value[1].Coeffs[0][0] = 1
	outputBefore := output.CopyNew()
	err := eval.ApplyEvaluationKey(badInput, evaluationKey, output)
	require.Error(t, err)
	require.Equal(t, outputBefore, output, "failed key application must not modify output")
	require.Equal(t, uint64(1), badInput.Value[1].Coeffs[0][0], "failed key application must not modify input")

	// Missing authoritative storage is not evidence of zero and must fail closed.
	malformed := input.CopyNew()
	malformed.Value[1].Coeffs[0] = nil
	outputBefore = output.CopyNew()
	err = eval.ApplyEvaluationKey(malformed, evaluationKey, output)
	require.Error(t, err)
	require.Equal(t, outputBefore, output, "malformed c1 backing must not modify the output")
	malformedKey := &rlwe.EvaluationKey{Layout: rlwe.KeyLayoutFast, GadgetCiphertext: rlwe.GadgetCiphertext{Layout: rlwe.KeyLayoutFast}}
	err = eval.ApplyEvaluationKey(input, malformedKey, output)
	require.Error(t, err)
	require.Equal(t, outputBefore, output, "malformed Fast key must not modify the output")

	level := input.Level()
	levelP := evaluationKey.LevelP()
	cx := params.RingQ().AtLevel(level).NewPoly()
	qpOut := rlwe.NewElementExtended(params, 1, level, levelP)
	require.NoError(t, eval.GadgetProductLazy(level, cx, &evaluationKey.GadgetCiphertext, qpOut))
	require.Equal(t, uint64(1), eval.fastQPCompat.counterSnapshot().gadgetProductElisions)
	for component := range qpOut.Value {
		require.True(t, fastQPolyIsZero(qpOut.Value[component].Q, level+1))
		require.True(t, fastQPolyIsZero(qpOut.Value[component].P, levelP+1))
	}

	qpBefore := qpOut.CopyNew()
	cx.Coeffs[0][0] = 1
	err = eval.GadgetProductLazy(level, cx, &evaluationKey.GadgetCiphertext, qpOut)
	require.Error(t, err)
	require.Equal(t, qpBefore, qpOut, "nonzero key-switch input must not modify QP output")
	cx.Zero()
	malformedGadget := &rlwe.GadgetCiphertext{Layout: rlwe.KeyLayoutFast}
	err = eval.GadgetProductLazy(level, cx, malformedGadget, qpOut)
	require.Error(t, err)
	require.Equal(t, qpBefore, qpOut, "malformed gadget storage must not modify QP output")

	beta := params.BaseRNSDecompositionVectorSize(level, levelP)
	decomposition := make([]ringqp.Poly, beta)
	for i := range decomposition {
		decomposition[i] = params.RingQP().AtLevel(level, levelP).NewPoly()
	}
	hoistedOut := rlwe.NewElementExtended(params, 1, level, levelP)
	require.NoError(t, eval.GadgetProductHoistedLazy(level, decomposition, &evaluationKey.GadgetCiphertext, hoistedOut))
	require.Equal(t, uint64(1), eval.fastQPCompat.counterSnapshot().hoistedProductElisions)
	for component := range hoistedOut.Value {
		require.True(t, fastQPolyIsZero(hoistedOut.Value[component].Q, level+1))
		require.True(t, fastQPolyIsZero(hoistedOut.Value[component].P, levelP+1))
	}

	legacyOut := rlwe.NewCiphertext(params, 1, level)
	eval.GadgetProductHoisted(level, decomposition, &evaluationKey.GadgetCiphertext, legacyOut)
	require.NoError(t, eval.fastQPCompat.peek(), "valid void GadgetProductHoisted must not defer an error")
	require.True(t, fastQPolyIsZero(legacyOut.Value[0], level+1))
	require.True(t, fastQPolyIsZero(legacyOut.Value[1], level+1))
	require.Equal(t, uint64(2), eval.fastQPCompat.counterSnapshot().hoistedProductElisions)
}

func TestFastQPCompatibilityFullQRescaleIsGatedAndMatchesRNSOracle(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	galEl := params.GaloisElementForRotation(1)
	gk := kgen.GenGaloisKeyNew(galEl, sk)
	eval := NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gk))
	encoder := NewEncoder(params)
	encryptor := rlwe.NewEncryptor(params, sk)
	input := fastRotateTestEncrypt(t, params, encoder, encryptor, fastRotateTestValues(4), 4, params.MaxLevel(), true)
	require.Equal(t, params.MaxLevel()+1, len(input.Value[0].Coeffs), "fixture must carry complete logical Q rows")

	// Before the verified high-level QP boundary, ordinary Fast Rescale stays
	// on the compact Q-prefix path and does not arm the full-Q fallback.
	compactInput := input.CopyNew()
	compactOut := rlwe.NewCiphertext(params, 1, compactInput.Level()-1)
	require.NoError(t, eval.Rescale(compactInput, compactOut))
	require.False(t, eval.fastQPCompat.isArmed())
	for row := fastcore.MaxQPrefixWidth; row <= compactOut.Level(); row++ {
		require.Empty(t, compactOut.Value[0].Coeffs[row])
	}

	beta := params.BaseRNSDecompositionVectorSize(input.Level(), gk.LevelP())
	decomposition := make([]ringqp.Poly, beta)
	for i := range decomposition {
		decomposition[i] = params.RingQP().AtLevel(input.Level(), gk.LevelP()).NewPoly()
	}
	// DecomposeNTT may inspect only the authoritative Q-prefix; dormant rows
	// are deliberately unallocated at this API boundary.
	compactC1 := input.CopyNew()
	for row := fastcore.MaxQPrefixWidth; row <= compactC1.Level(); row++ {
		compactC1.Value[1].Coeffs[row] = nil
	}
	eval.DecomposeNTT(input.Level(), gk.LevelP(), gk.LevelP()+1, compactC1.Value[1], input.IsNTT, decomposition)
	qp := rlwe.NewElementExtended(params, 1, input.Level(), gk.LevelP())
	require.NoError(t, eval.AutomorphismHoistedLazy(input.Level(), input, decomposition, galEl, qp))
	require.True(t, eval.fastQPCompat.isArmed())
	require.Equal(t, uint64(1), eval.fastQPCompat.counterSnapshot().decomposeElisions)
	modDown := rlwe.NewCiphertext(params, 1, input.Level())
	for component := range modDown.Value {
		eval.ModDownQPtoQNTT(input.Level(), gk.LevelP(), qp.Value[component].Q, qp.Value[component].P, modDown.Value[component])
	}
	require.NoError(t, eval.fastQPCompat.peek(), "verified high-Q ModDown must not defer an error")
	ringQFull := params.RingQ().AtLevel(input.Level())
	rotatedInput := ringQFull.NewPoly()
	ringQFull.AutomorphismNTTWithIndex(input.Value[0], eval.AutomorphismIndex(galEl), rotatedInput)
	for row := 0; row <= input.Level(); row++ {
		require.Equal(t, rotatedInput.Coeffs[row], modDown.Value[0].Coeffs[row], "full-Q ModDown row %d", row)
	}
	require.True(t, fastQPolyIsZero(modDown.Value[1], input.Level()+1))
	badP := qp.Value[0].P.CopyNew()
	badP.Coeffs[0][0] = 1
	modDownBefore := modDown.CopyNew()
	inputBeforeModDown := input.CopyNew()
	eval.ModDownQPtoQNTT(input.Level(), gk.LevelP(), qp.Value[0].Q, *badP, modDown.Value[0])
	require.Equal(t, modDownBefore, modDown, "rejected nonzero P input must leave the existing output unchanged")
	require.Error(t, eval.Trace(input, params.LogMaxSlots(), input), "void ModDown errors must surface at the next checked boundary")
	require.Equal(t, inputBeforeModDown, input, "failed ModDown must not mutate its source ciphertext")
	badQ := qp.Value[0].Q.CopyNew()
	q4 := params.RingQ().SubRings[fastcore.MaxQPrefixWidth].Modulus
	if badQ.Coeffs[fastcore.MaxQPrefixWidth][0]+1 == q4 {
		badQ.Coeffs[fastcore.MaxQPrefixWidth][0] = 0
	} else {
		badQ.Coeffs[fastcore.MaxQPrefixWidth][0]++
	}
	eval.ModDownQPtoQNTT(input.Level(), gk.LevelP(), *badQ, qp.Value[0].P, modDown.Value[0])
	require.Equal(t, modDownBefore, modDown, "incoherent dormant Q row must not reach full-RNS ModDown")
	require.Error(t, eval.Trace(input, params.LogMaxSlots(), input), "dormant-row mismatch must surface at the next checked boundary")

	got := rlwe.NewCiphertext(params, 1, input.Level()-1)
	inputBefore := input.CopyNew()
	require.NoError(t, eval.Rescale(input, got))
	require.Equal(t, uint64(1), eval.fastQPCompat.counterSnapshot().fullQRescales)
	require.Equal(t, inputBefore, input, "out-of-place compatibility Rescale must not modify input")
	require.Equal(t, input.Level()-1, got.Level())
	require.Equal(t, input.Scale.Div(rlwe.NewScale(params.Q()[input.Level()])), got.Scale)
	for component := range got.Value {
		for row := 0; row <= got.Level(); row++ {
			require.Len(t, got.Value[component].Coeffs[row], params.N(), "full-Q row %d component %d", row, component)
		}
	}
	for row := fastcore.MaxQPrefixWidth; row <= got.Level(); row++ {
		require.NotEmpty(t, got.Value[0].Coeffs[row], "compatibility Rescale must preserve non-prefix Q rows")
	}

	ringQ := params.RingQ().AtLevel(input.Level())
	oracle := rlwe.NewCiphertext(params, 1, input.Level()-1)
	buffer := ringQ.NewPoly()
	for component := range oracle.Value {
		ringQ.DivRoundByLastModulusManyNTT(1, input.Value[component], buffer, oracle.Value[component])
	}
	*oracle.MetaData = *input.MetaData
	oracle.Scale = input.Scale.Div(rlwe.NewScale(params.Q()[input.Level()]))
	for component := range oracle.Value {
		for row := 0; row <= oracle.Level(); row++ {
			require.Equal(t, oracle.Value[component].Coeffs[row], got.Value[component].Coeffs[row], "RNS oracle component %d row %d", component, row)
		}
	}
	require.True(t, got.Scale.Equal(oracle.Scale))
}
