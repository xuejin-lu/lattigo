package ckks

import (
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

type addSubCoreCall struct {
	rows int
	sub  bool
}

type recordingAddSubCore struct {
	delegate fastcore.AddSubCore
	calls    []addSubCoreCall
}

type t3ScalarCoreCall struct {
	rows           int
	levelsConsumed int
	real           string
	imag           string
	outputScale    string
}

type recordingT3ScalarCore struct {
	fastcore.ScalarCore
	calls []t3ScalarCoreCall
}

func (recorder *recordingT3ScalarCore) MulRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows, levelsConsumed int, scalar *bignum.Complex) error {
	err := recorder.ScalarCore.MulRows(ringQ, op0, opOut, rows, levelsConsumed, scalar)
	call := t3ScalarCoreCall{rows: rows, levelsConsumed: levelsConsumed, real: scalar[0].Text('f', 0), imag: scalar[1].Text('f', 0), outputScale: opOut.Scale.Value.Text('f', 0)}
	recorder.calls = append(recorder.calls, call)
	return err
}

type t3AddSubCoreCall struct {
	rows       int
	sub        bool
	op0Scale   string
	op1Scale   string
	opOutScale string
}

type recordingT3AddSubCore struct {
	fastcore.AddSubCore
	calls []t3AddSubCoreCall
}

func (recorder *recordingT3AddSubCore) ApplyRows(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, rows int, sub bool) error {
	call := t3AddSubCoreCall{rows: rows, sub: sub, op0Scale: op0.Scale.Value.Text('f', 0), op1Scale: op1.Scale.Value.Text('f', 0), opOutScale: opOut.Scale.Value.Text('f', 0)}
	recorder.calls = append(recorder.calls, call)
	return recorder.AddSubCore.ApplyRows(ringQ, op0, op1, opOut, rows, sub)
}

func newFastBTSScaleAlignmentCiphertext(params Parameters, degree, level int, scale rlwe.Scale, seed uint64) *rlwe.Ciphertext {
	ct := fastcore.NewCompactCiphertext(params, degree, level)
	ct.Scale = scale
	ct.IsBatched = true
	ct.LogDimensions = ring.Dimensions{Cols: 15}
	rows, _ := fastcore.QPrefixWidth(level)
	for row := 0; row < rows; row++ {
		modulus := params.RingQ().SubRings[row].Modulus
		for i := 0; i < 12; i++ {
			ct.Value[0].Coeffs[row][i] = (seed + uint64(row*17+i*3+1)) % modulus
		}
	}
	return ct
}

func requireFastBTSScaleAlignmentRowsUnchanged(t *testing.T, before, after *rlwe.Ciphertext, rows int) {
	t.Helper()
	require.Equal(t, before.Degree(), after.Degree())
	require.Equal(t, before.Level(), after.Level())
	for component := range before.Value {
		for row := 0; row < rows; row++ {
			require.Len(t, after.Value[component].Coeffs[row], len(before.Value[component].Coeffs[row]), "active row backing changed at component=%d q%d", component, row)
			for coefficient, value := range before.Value[component].Coeffs[row] {
				if after.Value[component].Coeffs[row][coefficient] != value {
					t.Fatalf("input changed at component=%d q%d coefficient=%d", component, row, coefficient)
				}
			}
		}
		for row := rows; row <= before.Level(); row++ {
			require.Len(t, after.Value[component].Coeffs[row], len(before.Value[component].Coeffs[row]), "dormant row backing changed at component=%d q%d", component, row)
		}
	}
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

func TestPublicFastZeroSecretT3DifferentScaleSubUsesQPrefixCores(t *testing.T) {
	params := btsCoverageLogN16Parameters(t)
	const level = 25
	q26 := params.Q()[26]
	lowScale := rlwe.NewScale(uint64(1) << 42)
	numerator := new(big.Int).Lsh(big.NewInt(1), 126)
	highScale := rlwe.NewScale(numerator).Div(rlwe.NewScale(q26))
	op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, highScale, 19)
	// BTS T3 multiplies T1(L26) by T2(L25), so the subtraction's low-scale
	// right operand remains at L26 while the in-place output is L25.
	op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level+1, lowScale, 71)
	op0Before := op0.CopyNew()
	op1Before := op1.CopyNew()

	evaluator := NewEvaluator(params, nil)
	scalarRecorder := &recordingT3ScalarCore{ScalarCore: fastcore.NewScalarWorkspace(params.RingQ(), fastcore.MaxQPrefixWidth)}
	addSubRecorder := &recordingT3AddSubCore{AddSubCore: fastcore.NewAddSubWorkspace()}
	evaluator.scalarCore = scalarRecorder
	evaluator.addSubCore = addSubRecorder

	// This proves modular semantics and dispatch for the observed experimental
	// path only. No trusted coefficient-bound certificate is asserted:
	// CAPACITY_UNPROVEN_EXPERIMENTAL.
	ratio := highScale.Div(lowScale).Value
	k, _ := ratio.Int(nil) // Standard's positive integer truncation semantics.
	require.True(t, k.Sign() > 0)
	exactRatio := new(big.Int).Lsh(big.NewInt(1), 84)
	exactRatio.Div(exactRatio, new(big.Int).SetUint64(q26))
	require.Equal(t, exactRatio, k, "the synthetic T3 metadata must retain the source-defined exact integer ratio")
	require.NoError(t, evaluator.Sub(op0, op1, op0))

	require.Len(t, scalarRecorder.calls, 1)
	require.Equal(t, fastcore.MaxQPrefixWidth, scalarRecorder.calls[0].rows)
	require.Equal(t, params.LevelsConsumedPerRescaling(), scalarRecorder.calls[0].levelsConsumed)
	require.Equal(t, k.String(), scalarRecorder.calls[0].real)
	require.Equal(t, "0", scalarRecorder.calls[0].imag)
	require.Equal(t, lowScale.Value.Text('f', 0), scalarRecorder.calls[0].outputScale, "scalar kernel must finish with the source Scale before metadata alignment")
	require.Equal(t, []t3AddSubCoreCall{{rows: fastcore.MaxQPrefixWidth, sub: true, op0Scale: highScale.Value.Text('f', 0), op1Scale: highScale.Value.Text('f', 0), opOutScale: highScale.Value.Text('f', 0)}}, addSubRecorder.calls)
	require.Equal(t, op1Before.Scale, op1.Scale, "right input metadata must remain unchanged")
	requireFastBTSScaleAlignmentRowsUnchanged(t, op1Before, op1, fastcore.MaxQPrefixWidth)

	require.Equal(t, level, op0.Level())
	require.Equal(t, level+1, op1.Level())
	require.Equal(t, 2, op0.Degree())
	require.True(t, highScale.Equal(op0.Scale))
	require.True(t, op0.IsNTT)
	require.False(t, op0.IsMontgomery)
	require.True(t, op0.IsBatched)
	require.Equal(t, ring.Dimensions{Cols: 15}, op0.LogDimensions)

	for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
		q := params.RingQ().SubRings[row].Modulus
		kModQ := new(big.Int).Mod(new(big.Int).Set(k), new(big.Int).SetUint64(q)).Uint64()
		for i, got := range op0.Value[0].Coeffs[row] {
			product := (kModQ * op1Before.Value[0].Coeffs[row][i]) % q
			want := (op0Before.Value[0].Coeffs[row][i] + q - product) % q
			require.Equal(t, want, got, "q%d coefficient %d", row, i)
		}
		for component := 1; component <= op0.Degree(); component++ {
			for _, coefficient := range op0.Value[component].Coeffs[row] {
				require.Zero(t, coefficient, "zero-secret component %d q%d", component, row)
			}
		}
	}
	for component := range op0.Value {
		for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
			require.Nil(t, op0.Value[component].Coeffs[row], "dormant q%d must remain unallocated", row)
		}
	}
}

func TestPublicFastZeroSecretT3RelinDegreeOneDifferentScaleSubUsesQPrefixCores(t *testing.T) {
	params := btsCoverageLogN16Parameters(t)
	const level = 25
	q26 := params.Q()[26]
	lowScale := rlwe.NewScale(uint64(1) << 42)
	highScale := rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 126)).Div(rlwe.NewScale(q26))
	op0 := newFastBTSScaleAlignmentCiphertext(params, 1, level, highScale, 23)
	op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level+1, lowScale, 79)
	op0Before, op1Before := op0.CopyNew(), op1.CopyNew()

	evaluator := NewEvaluator(params, nil)
	scalarRecorder := &recordingT3ScalarCore{ScalarCore: fastcore.NewScalarWorkspace(params.RingQ(), fastcore.MaxQPrefixWidth)}
	addSubRecorder := &recordingT3AddSubCore{AddSubCore: fastcore.NewAddSubWorkspace()}
	evaluator.scalarCore, evaluator.addSubCore = scalarRecorder, addSubRecorder

	ratio := highScale.Div(lowScale).Value
	k, _ := ratio.Int(nil)
	require.True(t, k.Sign() > 0)
	require.NoError(t, evaluator.Sub(op0, op1, op0))
	require.Len(t, scalarRecorder.calls, 1)
	require.Equal(t, fastcore.MaxQPrefixWidth, scalarRecorder.calls[0].rows)
	require.Equal(t, k.String(), scalarRecorder.calls[0].real)
	require.Equal(t, "0", scalarRecorder.calls[0].imag)
	require.Equal(t, []t3AddSubCoreCall{{rows: fastcore.MaxQPrefixWidth, sub: true, op0Scale: highScale.Value.Text('f', 0), op1Scale: highScale.Value.Text('f', 0), opOutScale: highScale.Value.Text('f', 0)}}, addSubRecorder.calls)
	require.Equal(t, level, op0.Level())
	require.Equal(t, level+1, op1.Level())
	require.Equal(t, 1, op0.Degree())
	require.True(t, highScale.Equal(op0.Scale))
	require.True(t, op0.IsNTT)
	require.False(t, op0.IsMontgomery)
	require.True(t, op0.IsBatched)
	require.Equal(t, ring.Dimensions{Cols: 15}, op0.LogDimensions)
	requireFastBTSScaleAlignmentRowsUnchanged(t, op1Before, op1, fastcore.MaxQPrefixWidth)

	for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
		q := params.RingQ().SubRings[row].Modulus
		kModQ := new(big.Int).Mod(new(big.Int).Set(k), new(big.Int).SetUint64(q)).Uint64()
		for i, got := range op0.Value[0].Coeffs[row] {
			product := (kModQ * op1Before.Value[0].Coeffs[row][i]) % q
			want := (op0Before.Value[0].Coeffs[row][i] + q - product) % q
			require.Equal(t, want, got, "q%d coefficient %d", row, i)
		}
		for _, coefficient := range op0.Value[1].Coeffs[row] {
			require.Zero(t, coefficient, "zero-secret c1 q%d", row)
		}
	}
	for _, ct := range []*rlwe.Ciphertext{op0, op1} {
		for component := range ct.Value {
			for row := fastcore.MaxQPrefixWidth; row <= ct.Level(); row++ {
				require.Nil(t, ct.Value[component].Coeffs[row], "dormant q%d must remain unallocated", row)
			}
		}
	}
}

func TestPublicFastZeroSecretT3DifferentScaleSubRejectsMalformedOrUnsupportedInputs(t *testing.T) {
	params := btsCoverageLogN16Parameters(t)
	const level = 25
	low := rlwe.NewScale(uint64(1) << 42)
	q26 := params.Q()[26]
	high := rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 126)).Div(rlwe.NewScale(q26))
	evaluator := NewEvaluator(params, nil)
	scalarRecorder := &recordingT3ScalarCore{ScalarCore: fastcore.NewScalarWorkspace(params.RingQ(), fastcore.MaxQPrefixWidth)}
	addSubRecorder := &recordingT3AddSubCore{AddSubCore: fastcore.NewAddSubWorkspace()}
	evaluator.scalarCore = scalarRecorder
	evaluator.addSubCore = addSubRecorder

	t.Run("nonzero c1", func(t *testing.T) {
		op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, high, 5)
		op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level, low, 9)
		op1.Value[1].Coeffs[0][0] = 1
		before0, before1 := op0.CopyNew(), op1.CopyNew()
		err := evaluator.Sub(op0, op1, op0)
		require.Error(t, err, "compact nonzero-c1 input must not reach the full-Q evaluator")
		requireFastBTSScaleAlignmentRowsUnchanged(t, before0, op0, fastcore.MaxQPrefixWidth)
		require.Equal(t, before0.Scale, op0.Scale)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before1, op1, fastcore.MaxQPrefixWidth)
	})

	t.Run("nonzero c2", func(t *testing.T) {
		op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, high, 7)
		op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level, low, 11)
		op0.Value[2].Coeffs[2][3] = 1
		before0, before1 := op0.CopyNew(), op1.CopyNew()
		err := evaluator.Sub(op0, op1, op0)
		require.Error(t, err, "compact nonzero-c2 input must not reach the full-Q evaluator")
		requireFastBTSScaleAlignmentRowsUnchanged(t, before0, op0, fastcore.MaxQPrefixWidth)
		require.Equal(t, before0.Scale, op0.Scale)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before1, op1, fastcore.MaxQPrefixWidth)
	})

	t.Run("missing authoritative q3 backing", func(t *testing.T) {
		op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, high, 13)
		op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level, low, 17)
		op1.Value[0].Coeffs[3] = nil
		before0, before1 := op0.CopyNew(), op1.CopyNew()
		err := evaluator.Sub(op0, op1, op0)
		require.Error(t, err)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before0, op0, fastcore.MaxQPrefixWidth)
		require.Equal(t, before0.Scale, op0.Scale)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before1, op1, fastcore.MaxQPrefixWidth)
	})

	t.Run("noncanonical residue", func(t *testing.T) {
		op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, high, 47)
		op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level, low, 53)
		op1.Value[0].Coeffs[0][0] = params.RingQ().SubRings[0].Modulus
		before0, before1 := op0.CopyNew(), op1.CopyNew()
		err := evaluator.Sub(op0, op1, op0)
		require.ErrorContains(t, err, "not a canonical residue")
		requireFastBTSScaleAlignmentRowsUnchanged(t, before0, op0, fastcore.MaxQPrefixWidth)
		require.Equal(t, before0.Scale, op0.Scale)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before1, op1, fastcore.MaxQPrefixWidth)
	})

	t.Run("lower scale on left is not aligned", func(t *testing.T) {
		op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, low, 23)
		op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level, high, 29)
		before0, before1 := op0.CopyNew(), op1.CopyNew()
		err := evaluator.Sub(op0, op1, op0)
		require.Error(t, err, "unsupported compact scale direction must fail before mutation")
		requireFastBTSScaleAlignmentRowsUnchanged(t, before0, op0, fastcore.MaxQPrefixWidth)
		require.Equal(t, before0.Scale, op0.Scale)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before1, op1, fastcore.MaxQPrefixWidth)
	})

	t.Run("separate output alias is not aligned", func(t *testing.T) {
		op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, high, 31)
		op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level, low, 37)
		opOut := fastcore.NewCompactCiphertext(params, 2, level)
		opOut.Scale = high
		before0, before1, beforeOut := op0.CopyNew(), op1.CopyNew(), opOut.CopyNew()
		err := evaluator.Sub(op0, op1, opOut)
		require.Error(t, err, "unsupported alias direction must fail before mutation")
		requireFastBTSScaleAlignmentRowsUnchanged(t, before0, op0, fastcore.MaxQPrefixWidth)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before1, op1, fastcore.MaxQPrefixWidth)
		requireFastBTSScaleAlignmentRowsUnchanged(t, beforeOut, opOut, fastcore.MaxQPrefixWidth)
		require.Equal(t, beforeOut.Scale, opOut.Scale)
	})

	t.Run("unequal logical Levels are not aligned", func(t *testing.T) {
		op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, high, 41)
		op1 := newFastBTSScaleAlignmentCiphertext(params, 1, level-1, low, 43)
		before0, before1 := op0.CopyNew(), op1.CopyNew()
		err := evaluator.Sub(op0, op1, op0)
		require.Error(t, err, "unsupported logical Levels must fail before mutation")
		requireFastBTSScaleAlignmentRowsUnchanged(t, before0, op0, fastcore.MaxQPrefixWidth)
		require.Equal(t, before0.Scale, op0.Scale)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before1, op1, fastcore.MaxQPrefixWidth)
	})

	t.Run("unsupported right Degree is not aligned", func(t *testing.T) {
		op0 := newFastBTSScaleAlignmentCiphertext(params, 2, level, high, 59)
		op1 := newFastBTSScaleAlignmentCiphertext(params, 0, level, low, 61)
		before0, before1 := op0.CopyNew(), op1.CopyNew()
		err := evaluator.Sub(op0, op1, op0)
		require.Error(t, err, "unsupported component Degree must fail before mutation")
		requireFastBTSScaleAlignmentRowsUnchanged(t, before0, op0, fastcore.MaxQPrefixWidth)
		require.Equal(t, before0.Scale, op0.Scale)
		requireFastBTSScaleAlignmentRowsUnchanged(t, before1, op1, fastcore.MaxQPrefixWidth)
	})

	require.Empty(t, scalarRecorder.calls, "rejected inputs must not reach the scalar kernel")
	require.Empty(t, addSubRecorder.calls, "rejected inputs must not reach the Add/Sub kernel")
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

	// Public MulRelin consumes exactly q0..q3 and keeps logical Level 5;
	// dormant q4/q5 stay absent rather than being materialized or read.
	mulOut, err := evaluator.MulRelinNew(sum, difference)
	require.NoError(t, err)
	require.Equal(t, level, mulOut.Level())
	require.Equal(t, 1, mulOut.Degree())
	require.True(t, mulOut.Scale.Equal(sum.Scale.Mul(difference.Scale)))
	for component := range mulOut.Value {
		for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
			require.Len(t, mulOut.Value[component].Coeffs[row], params.N())
		}
		for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
			require.Nil(t, mulOut.Value[component].Coeffs[row], "dormant q%d must remain unmaterialized", row)
		}
	}
	for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
		for _, coefficient := range mulOut.Value[1].Coeffs[row] {
			require.Zero(t, coefficient)
		}
	}

	// The same scratch-first core remains safe when output aliases either input.
	mulAlias := sum.CopyNew()
	require.NoError(t, evaluator.MulRelin(mulAlias, difference, mulAlias))
	require.Equal(t, level, mulAlias.Level())
	require.Equal(t, 1, mulAlias.Degree())
	for row := fastcore.MaxQPrefixWidth; row <= level; row++ {
		require.Nil(t, mulAlias.Value[0].Coeffs[row])
		require.Nil(t, mulAlias.Value[1].Coeffs[row])
	}

	// Level 5 remains structural-only: q4/q5 stay dormant through the
	// existing Rescale and Rotate consumers, with no numerical claim here.
	rescaled := NewCiphertext(params, 1, level-1)
	require.NoError(t, evaluator.Rescale(mulOut, rescaled))
	require.Equal(t, level-1, rescaled.Level())
	require.True(t, mulOut.Scale.Div(rlwe.NewScale(params.Q()[level])).Equal(rescaled.Scale))
	for component := range rescaled.Value {
		for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
			require.Len(t, rescaled.Value[component].Coeffs[row], params.N())
		}
		for row := fastcore.MaxQPrefixWidth; row <= rescaled.Level(); row++ {
			require.Nil(t, rescaled.Value[component].Coeffs[row], "dormant q%d must remain unmaterialized", row)
		}
	}
	for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
		for _, coefficient := range rescaled.Value[1].Coeffs[row] {
			require.Zero(t, coefficient)
		}
	}

	rotated, err := evaluator.RotateNew(rescaled, 1)
	require.NoError(t, err)
	require.Equal(t, level-1, rotated.Level())
	for component := range rotated.Value {
		for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
			require.Len(t, rotated.Value[component].Coeffs[row], params.N())
		}
		for row := fastcore.MaxQPrefixWidth; row <= rotated.Level(); row++ {
			require.Nil(t, rotated.Value[component].Coeffs[row], "dormant q%d must remain unmaterialized", row)
		}
	}
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
