package ckks

import (
	"math"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

type addSubCoreCall struct {
	rows int
	sub  bool
}

type recordingAddSubCore struct {
	delegate fastcore.AddSubCore
	calls    []addSubCoreCall
}

func (recorder *recordingAddSubCore) ApplyRows(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, rows int, sub bool) error {
	recorder.calls = append(recorder.calls, addSubCoreCall{rows: rows, sub: sub})
	return recorder.delegate.ApplyRows(ringQ, op0, op1, opOut, rows, sub)
}

func fastAddSubTestValues(seed int) []complex128 {
	values := make([]complex128, 1<<fastRotateTestLogSlots)
	for i := range values {
		values[i] = complex(float64((i+seed)%13-6)/32, float64((3*i+seed)%9-4)/64)
	}
	return values
}

func encryptFastAddSubTestValues(t *testing.T, params Parameters, encoder *Encoder, encryptor *rlwe.Encryptor, values []complex128, level int) *rlwe.Ciphertext {
	t.Helper()
	pt := NewPlaintext(params, level)
	pt.Scale = rlwe.NewScale(math.Exp2(fastRotateTestLogScale))
	pt.LogDimensions = ring.Dimensions{Cols: fastRotateTestLogSlots}
	require.NoError(t, encoder.Encode(values, pt))
	ct, err := encryptor.EncryptNew(pt)
	require.NoError(t, err)
	return ct
}

func decodeFastAddSubTestValues(t *testing.T, params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, ct *rlwe.Ciphertext) []complex128 {
	t.Helper()
	values := make([]complex128, 1<<fastRotateTestLogSlots)
	require.NoError(t, encoder.Decode(decryptor.DecryptNew(ct), values))
	return values
}

func requireComplexValuesNear(t *testing.T, got, want []complex128, tolerance float64) {
	t.Helper()
	require.Len(t, got, len(want))
	for i := range want {
		require.LessOrEqual(t, math.Abs(real(got[i]-want[i])), tolerance, "real slot %d", i)
		require.LessOrEqual(t, math.Abs(imag(got[i]-want[i])), tolerance, "imag slot %d", i)
	}
}

func TestPublicFastZeroSecretAddSubUsesSharedQPrefixCoreAndDecodes(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	decryptor := rlwe.NewDecryptor(params, sk)
	values0, values1 := fastAddSubTestValues(1), fastAddSubTestValues(4)
	evaluator := NewEvaluator(params, nil)
	recorder := &recordingAddSubCore{delegate: fastcore.NewAddSubWorkspace()}
	evaluator.addSubCore = recorder

	for _, level := range []int{1, 3} {
		ct0 := encryptFastAddSubTestValues(t, params, encoder, encryptor, values0, level)
		ct1 := encryptFastAddSubTestValues(t, params, encoder, encryptor, values1, level)
		sum, err := evaluator.AddNew(ct0, ct1)
		require.NoError(t, err)
		require.Equal(t, level, sum.Level())
		require.Equal(t, 1, sum.Degree())
		require.True(t, ct0.Scale.Equal(sum.Scale))
		require.Equal(t, level+1, len(sum.Value[0].Coeffs))
		for row := 0; row <= level; row++ {
			require.Len(t, sum.Value[0].Coeffs[row], params.N())
			require.Len(t, sum.Value[1].Coeffs[row], params.N())
			for _, coefficient := range sum.Value[1].Coeffs[row] {
				require.Zero(t, coefficient)
			}
		}
		wantSum := make([]complex128, len(values0))
		for i := range wantSum {
			wantSum[i] = values0[i] + values1[i]
		}
		requireComplexValuesNear(t, decodeFastAddSubTestValues(t, params, encoder, decryptor, sum), wantSum, 1e-6)
		require.Equal(t, addSubCoreCall{rows: min(level+1, fastcore.MaxQPrefixWidth)}, recorder.calls[len(recorder.calls)-1])

		// Public Sub consumes the compact Add output through the same core.
		difference, err := evaluator.SubNew(sum, ct1)
		require.NoError(t, err)
		require.Equal(t, level, difference.Level())
		require.True(t, sum.Scale.Equal(difference.Scale))
		requireComplexValuesNear(t, decodeFastAddSubTestValues(t, params, encoder, decryptor, difference), values0, 1e-6)
		require.Equal(t, addSubCoreCall{rows: min(level+1, fastcore.MaxQPrefixWidth), sub: true}, recorder.calls[len(recorder.calls)-1])
	}
	require.Len(t, recorder.calls, 4)
}

func TestPublicFastZeroSecretAddSubCompactsAbovePrefixAndRejectsUnsafeConsumers(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	evaluator := NewEvaluator(params, nil)
	recorder := &recordingAddSubCore{delegate: fastcore.NewAddSubWorkspace()}
	evaluator.addSubCore = recorder
	const level = 5
	values0, values1 := fastAddSubTestValues(2), fastAddSubTestValues(7)
	ct0 := encryptFastAddSubTestValues(t, params, encoder, encryptor, values0, level)
	ct1 := encryptFastAddSubTestValues(t, params, encoder, encryptor, values1, level)

	// Poison rows above the Q-prefix. They are not authoritative or read by the
	// compact public bridge.
	for _, ct := range []*rlwe.Ciphertext{ct0, ct1} {
		for component := range ct.Value {
			for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
				for i := range ct.Value[component].Coeffs[row] {
					ct.Value[component].Coeffs[row][i] = ^uint64(0)
				}
			}
		}
	}

	sum, err := evaluator.AddNew(ct0, ct1)
	require.NoError(t, err)
	require.Equal(t, level, sum.Level(), "logical CKKS Level remains unchanged")
	for component := range sum.Value {
		require.Len(t, sum.Value[component].Coeffs, level+1)
		for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
			require.Len(t, sum.Value[component].Coeffs[row], params.N())
		}
		for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
			require.Nil(t, sum.Value[component].Coeffs[row], "dormant q%d must remain unmaterialized", row)
		}
	}
	for component := 0; component <= 1; component++ {
		for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
			want := make([]uint64, params.N())
			params.RingQ().SubRings[row].Add(ct0.Value[component].Coeffs[row], ct1.Value[component].Coeffs[row], want)
			require.Equal(t, want, sum.Value[component].Coeffs[row], "component=%d q%d", component, row)
		}
	}
	require.Equal(t, addSubCoreCall{rows: fastcore.MaxQPrefixWidth}, recorder.calls[0])

	difference, err := evaluator.SubNew(sum, ct1)
	require.NoError(t, err)
	for component := range difference.Value {
		for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
			require.Equal(t, ct0.Value[component].Coeffs[row], difference.Value[component].Coeffs[row])
		}
		for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
			require.Nil(t, difference.Value[component].Coeffs[row])
		}
	}
	require.Equal(t, addSubCoreCall{rows: fastcore.MaxQPrefixWidth, sub: true}, recorder.calls[1])

	// Existing full-active-Q public consumers must reject this compact result
	// before reading q4 or mutating their outputs.
	mulOut := NewCiphertext(params, 1, level)
	mulOutBefore := mulOut.CopyNew()
	err = evaluator.MulRelin(sum, difference, mulOut)
	require.ErrorContains(t, err, "q4")
	require.Equal(t, mulOutBefore.Value, mulOut.Value)
	require.Equal(t, mulOutBefore.MetaData, mulOut.MetaData)

	rescaleOut := NewCiphertext(params, 1, level-1)
	rescaleOutBefore := rescaleOut.CopyNew()
	err = evaluator.Rescale(sum, rescaleOut)
	require.ErrorContains(t, err, "q4")
	require.Equal(t, rescaleOutBefore.Value, rescaleOut.Value)
	require.Equal(t, rescaleOutBefore.MetaData, rescaleOut.MetaData)

	_, err = evaluator.RotateNew(sum, 1)
	require.ErrorContains(t, err, "q4")
	require.Equal(t, 2, len(recorder.calls), "unsupported consumers must not fall back through the Add/Sub core")
}

func TestPublicFastZeroSecretAddFallsBackForFullCiphertextAndRejectsCompactGenericOperand(t *testing.T) {
	params := fastRotateTestParameters(t, ring.Standard)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	evaluator := NewEvaluator(params, nil)
	recorder := &recordingAddSubCore{delegate: fastcore.NewAddSubWorkspace()}
	evaluator.addSubCore = recorder
	level := 3
	values := fastAddSubTestValues(3)
	ct0 := encryptFastAddSubTestValues(t, params, encoder, encryptor, values, level)
	ct1 := encryptFastAddSubTestValues(t, params, encoder, encryptor, values, level)
	bad := ct0.CopyNew()
	bad.Value[1].Coeffs[0][0] = 1
	out := NewCiphertext(params, 1, level)
	err := evaluator.Add(bad, ct1, out)
	require.NoError(t, err, "full-Q ciphertexts outside the Fast zero-secret contract retain the generic CKKS path")
	require.Empty(t, recorder.calls, "invalid input must not dispatch")
	require.Equal(t, uint64(1), out.Value[1].Coeffs[0][0], "generic fallback must preserve the nonzero c1 contribution")

	badResidue := ct0.CopyNew()
	badResidue.Value[0].Coeffs[0][0] = params.Q()[0]
	residueOut := NewCiphertext(params, 1, level)
	residueOutBefore := residueOut.CopyNew()
	err = evaluator.Add(badResidue, ct1, residueOut)
	require.ErrorContains(t, err, "not a canonical residue")
	require.Equal(t, residueOutBefore.Value, residueOut.Value, "invalid specialized input must not mutate the output")
	require.Empty(t, recorder.calls, "invalid residue must not enter the Fast core")

	// The four-row cap is structural, not inferred from full backing: once a
	// public Fast result is compact, scalar Add must not enter generic full-Q.
	paramsLong := fastRotateTestParameters(t, ring.Standard)
	ctLong0 := encryptFastAddSubTestValues(t, paramsLong, NewEncoder(paramsLong), rlwe.NewEncryptor(paramsLong, rlwe.NewKeyGenerator(paramsLong).GenSecretKeyNew()), values, 5)
	ctLong1 := ctLong0.CopyNew()
	evaluatorLong := NewEvaluator(paramsLong, nil)
	recorderLong := &recordingAddSubCore{delegate: fastcore.NewAddSubWorkspace()}
	evaluatorLong.addSubCore = recorderLong
	compact, err := evaluatorLong.AddNew(ctLong0, ctLong1)
	require.NoError(t, err)
	compactOut := NewCiphertext(paramsLong, 1, 5)
	compactOutBefore := compactOut.CopyNew()
	err = evaluatorLong.Add(compact, int64(3), compactOut)
	require.ErrorContains(t, err, "missing active row q4")
	require.Equal(t, compactOutBefore.Value, compactOut.Value)
}
