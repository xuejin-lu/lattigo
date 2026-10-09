package ckks

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

type mulCoreCall struct {
	rows  int
	relin bool
}

type recordingMulCore struct {
	delegate fastcore.MulCore
	calls    []mulCoreCall
}

func (recorder *recordingMulCore) ApplyRows(ringQ *ring.Ring, op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], opOut *rlwe.Ciphertext, rows int, relin bool) error {
	recorder.calls = append(recorder.calls, mulCoreCall{rows: rows, relin: relin})
	return recorder.delegate.ApplyRows(ringQ, op0, op1, opOut, rows, relin)
}

func TestPublicFastZeroSecretMulRelinUsesSharedCoreAndDecodes(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	decryptor := rlwe.NewDecryptor(params, sk)
	evaluator := NewEvaluator(params, nil)
	recorder := &recordingMulCore{delegate: fastcore.NewMulWorkspace(params.RingQ(), fastcore.MaxQPrefixWidth)}
	evaluator.mulCore = recorder
	values0, values1 := fastAddSubTestValues(2), fastAddSubTestValues(6)

	for _, level := range []int{1, 3} {
		a := encryptFastAddSubTestValues(t, params, encoder, encryptor, values0, level)
		b := encryptFastAddSubTestValues(t, params, encoder, encryptor, values1, level)
		got, err := evaluator.MulRelinNew(a, b)
		require.NoError(t, err)
		require.Equal(t, level, got.Level())
		require.Equal(t, 1, got.Degree())
		require.True(t, got.Scale.Equal(a.Scale.Mul(b.Scale)))
		require.True(t, got.IsNTT)
		require.False(t, got.IsMontgomery)
		require.Equal(t, level+1, len(got.Value[0].Coeffs))
		for q := 0; q <= level; q++ {
			for _, coefficient := range got.Value[1].Coeffs[q] {
				require.Zero(t, coefficient)
			}
		}
		want := make([]complex128, len(values0))
		for i := range want {
			want[i] = values0[i] * values1[i]
		}
		requireComplexValuesNear(t, decodeFastAddSubTestValues(t, params, encoder, decryptor, got), want, 1e-4)
		require.Equal(t, mulCoreCall{rows: min(level+1, fastcore.MaxQPrefixWidth), relin: true}, recorder.calls[len(recorder.calls)-1])
	}
	require.Len(t, recorder.calls, 2)

	// A full-backed input with nonzero c1 must fail closed, before the shared
	// product core or any evaluation-key path can run.
	a := encryptFastAddSubTestValues(t, params, encoder, encryptor, values0, 3)
	b := encryptFastAddSubTestValues(t, params, encoder, encryptor, values1, 3)
	a.Value[1].Coeffs[0][0] = 1
	out := fastcore.NewCompactCiphertext(params, 1, 3)
	before, _ := out.MarshalBinary()
	require.ErrorContains(t, evaluator.MulRelin(a, b, out), "nonzero c1")
	after, _ := out.MarshalBinary()
	require.Equal(t, before, after)
	require.Len(t, recorder.calls, 2)
}

func TestPublicFastZeroSecretMulRelinNewRejectsMalformedCiphertextsWithoutPanic(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	evaluator := NewEvaluator(params, nil)
	valid := encryptFastAddSubTestValues(t, params, encoder, encryptor, fastAddSubTestValues(2), 3)
	other := encryptFastAddSubTestValues(t, params, encoder, encryptor, fastAddSubTestValues(6), 3)

	for _, malformed := range []string{"no_components", "one_component", "nil_metadata", "empty_poly", "empty_q0", "short_component_rows"} {
		t.Run(malformed, func(t *testing.T) {
			input, wantInput := valid.CopyNew(), valid.CopyNew()
			corruptFastMulRelinCiphertext(input, malformed)
			corruptFastMulRelinCiphertext(wantInput, malformed)

			got, err := evaluator.MulRelinNew(input, other)
			require.Error(t, err)
			require.Nil(t, got)
			require.Equal(t, wantInput.Value, input.Value)
			require.Equal(t, wantInput.MetaData, input.MetaData)
		})
	}
}

func TestPublicFastZeroSecretMulRelinMalformedInputLeavesInputsAndOutputUnchanged(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	evaluator := NewEvaluator(params, nil)
	op0 := encryptFastAddSubTestValues(t, params, encoder, encryptor, fastAddSubTestValues(2), 3)
	op1 := encryptFastAddSubTestValues(t, params, encoder, encryptor, fastAddSubTestValues(6), 3)
	want0, want1 := op0.CopyNew(), op1.CopyNew()
	corruptFastMulRelinCiphertext(op1, "empty_q0")
	corruptFastMulRelinCiphertext(want1, "empty_q0")

	output := fastcore.NewCompactCiphertext(params, 1, 3)
	output.Scale = rlwe.NewScale(1 << 27)
	output.MetaData.IsBatched = true
	for component := range output.Value {
		for row := 0; row < 4; row++ {
			for i := range output.Value[component].Coeffs[row] {
				output.Value[component].Coeffs[row][i] = uint64(19 + component + row + i)
			}
		}
	}
	wantOutput := output.CopyNew()

	err := evaluator.MulRelin(op0, op1, output)
	require.Error(t, err)
	require.Equal(t, want0.Value, op0.Value)
	require.Equal(t, want0.MetaData, op0.MetaData)
	require.Equal(t, want1.Value, op1.Value)
	require.Equal(t, want1.MetaData, op1.MetaData)
	require.Equal(t, wantOutput.Value, output.Value)
	require.Equal(t, wantOutput.MetaData, output.MetaData)
	require.True(t, wantOutput.Scale.Equal(output.Scale))
	require.Equal(t, wantOutput.IsNTT, output.IsNTT)
}

func corruptFastMulRelinCiphertext(ct *rlwe.Ciphertext, kind string) {
	switch kind {
	case "no_components":
		ct.Value = nil
	case "one_component":
		ct.Value = ct.Value[:1]
	case "nil_metadata":
		ct.MetaData = nil
	case "empty_poly":
		ct.Value[0].Coeffs = nil
	case "empty_q0":
		ct.Value[0].Coeffs[0] = nil
	case "short_component_rows":
		ct.Value[1].Coeffs = ct.Value[1].Coeffs[:2]
	default:
		panic("unknown malformed ciphertext case")
	}
}
