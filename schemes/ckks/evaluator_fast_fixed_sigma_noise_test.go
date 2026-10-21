package ckks

import (
	"encoding/hex"
	"math"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func fixedSigmaTestRescaledOutput(t *testing.T) (Parameters, *Evaluator, *rlwe.Ciphertext) {
	t.Helper()
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39, 40, 39, 40, 60, 60}, LogP: []int{61},
		LogDefaultScale: 45, RingType: ring.Standard,
	})
	require.NoError(t, err)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	evaluator := NewEvaluator(params, nil)
	x := fixedSigmaTestEncrypt(t, params, encoder, encryptor, fastAddSubTestValues(2), rlwe.NewScale(math.Exp2(45)))
	y := fixedSigmaTestEncrypt(t, params, encoder, encryptor, fastAddSubTestValues(6), rlwe.NewScale(params.Q()[5]))
	product, err := evaluator.MulRelinNew(x, y)
	require.NoError(t, err)
	output := NewCiphertext(params, 1, 4)
	require.NoError(t, evaluator.Rescale(product, output))
	require.Equal(t, 4, output.Level())
	require.True(t, output.IsNTT)
	require.False(t, output.IsMontgomery)
	require.True(t, output.Scale.Equal(rlwe.NewScale(math.Exp2(45))))
	return params, evaluator, output
}

func fixedSigmaTestEncrypt(t *testing.T, params Parameters, encoder *Encoder, encryptor *rlwe.Encryptor, values []complex128, scale rlwe.Scale) *rlwe.Ciphertext {
	t.Helper()
	plain := NewPlaintext(params, 5)
	plain.Scale, plain.IsNTT = scale, true
	plain.LogDimensions = ring.Dimensions{Cols: 4}
	require.NoError(t, encoder.Encode(values, plain))
	ct, err := encryptor.EncryptNew(plain)
	require.NoError(t, err)
	return ct
}

func TestFixedSigmaCoefficientNoiseZeroIsStrictNoOp(t *testing.T) {
	_, evaluator, output := fixedSigmaTestRescaledOutput(t)
	before := snapshotCiphertextValue(output)
	evidence, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(output, 0, 17)
	require.NoError(t, err)
	require.True(t, evidence.NoOp)
	require.Zero(t, evidence.SampleCount)
	require.Empty(t, evidence.CoefficientSHA256)
	require.Equal(t, before, snapshotCiphertextValue(output))

	// The zero-sigma path deliberately does not inspect the receiver or target.
	evidence, err = (Evaluator{}).InjectFixedSigmaCoefficientNoiseDiagnostic(nil, 0, 99)
	require.NoError(t, err)
	require.True(t, evidence.NoOp)
}

func TestFixedSigmaCoefficientNoiseIsDeterministicAndMapsOnePolynomialToEveryRow(t *testing.T) {
	params, evaluator, output := fixedSigmaTestRescaledOutput(t)
	const sigma, seed = 3.75, int64(20261022)
	first, second := output.CopyNew(), output.CopyNew()
	before := first.CopyNew()
	c1Before := snapshotCiphertextValueComponent(first, 1)
	dormantBefore := rowBackingLengths(first.Value[0].Coeffs[4])
	evidenceFirst, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(first, sigma, seed)
	require.NoError(t, err)
	evidenceSecond, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(second, sigma, seed)
	require.NoError(t, err)
	require.Equal(t, evidenceFirst, evidenceSecond)
	require.Equal(t, first.Value[0], second.Value[0])
	require.Equal(t, evidenceFirst.SampleCount, params.N())
	require.Equal(t, 4, evidenceFirst.Capacity.Rows)
	require.Equal(t, params.N(), evidenceFirst.Capacity.CheckedCoefficients)
	require.True(t, evidenceFirst.Capacity.StrictFit)
	require.NotEmpty(t, evidenceFirst.Capacity.PrefixProduct)
	require.NotEmpty(t, evidenceFirst.CoefficientSHA256)
	decodedHash, err := hex.DecodeString(evidenceFirst.CoefficientSHA256)
	require.NoError(t, err)
	require.Len(t, decodedHash, 32)

	// Only c0 changes. Ciphertext metadata, c1, and the dormant q4 backing stay
	// identical; the test never inspects dormant row contents.
	require.Equal(t, before.Degree(), first.Degree())
	require.Equal(t, before.Level(), first.Level())
	require.Equal(t, *before.MetaData, *first.MetaData)
	require.Equal(t, before.IsNTT, first.IsNTT)
	require.Equal(t, before.IsMontgomery, first.IsMontgomery)
	require.Equal(t, c1Before, snapshotCiphertextValueComponent(first, 1))
	require.Equal(t, dormantBefore, rowBackingLengths(first.Value[0].Coeffs[4]))
	require.Empty(t, dormantBefore)

	// Recreate the specified rounded Gaussian sequence, then verify the
	// inverse-transformed delta in every maintained residue row equals that
	// same signed integer polynomial modulo the row's q_i.
	rng := rand.New(rand.NewSource(seed))
	expected := make([]int64, params.N())
	for coefficient := range expected {
		expected[coefficient] = int64(math.Round(rng.NormFloat64() * sigma))
	}
	for row := 0; row < 4; row++ {
		subring := params.RingQ().SubRings[row]
		q := subring.Modulus
		deltaNTT := make([]uint64, params.N())
		for coefficient := range deltaNTT {
			after, prior := first.Value[0].Coeffs[row][coefficient], before.Value[0].Coeffs[row][coefficient]
			if after >= prior {
				deltaNTT[coefficient] = after - prior
			} else {
				deltaNTT[coefficient] = q - (prior - after)
			}
		}
		coefficientDelta := make([]uint64, params.N())
		subring.INTT(deltaNTT, coefficientDelta)
		for coefficient, got := range coefficientDelta {
			require.Equal(t, signedCoefficientResidue(expected[coefficient], q), got,
				"q%d coefficient %d must encode the same sampled signed integer", row, coefficient)
		}
	}
}

func TestFixedSigmaCoefficientNoiseRejectsNonzeroC1WithoutMutation(t *testing.T) {
	_, evaluator, output := fixedSigmaTestRescaledOutput(t)
	output.Value[1].Coeffs[1][0] = 1
	before := snapshotCiphertextValue(output)
	metadata := *output.MetaData
	_, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(output, 3.75, 9)
	require.ErrorContains(t, err, "c1=0")
	require.Equal(t, before, snapshotCiphertextValue(output))
	require.Equal(t, metadata, *output.MetaData)
}

func snapshotCiphertextValue(ct *rlwe.Ciphertext) [][][]uint64 {
	values := make([][][]uint64, len(ct.Value))
	for component := range ct.Value {
		values[component] = snapshotCiphertextValueComponent(ct, component)
	}
	return values
}

func snapshotCiphertextValueComponent(ct *rlwe.Ciphertext, component int) [][]uint64 {
	poly := ct.Value[component]
	if poly.Coeffs == nil {
		return nil
	}
	rows := make([][]uint64, len(poly.Coeffs))
	for row, values := range poly.Coeffs {
		if values != nil {
			rows[row] = append([]uint64{}, values...)
		}
	}
	return rows
}

func rowBackingLengths(values []uint64) []int {
	if len(values) == 0 {
		return nil
	}
	return []int{len(values)}
}
