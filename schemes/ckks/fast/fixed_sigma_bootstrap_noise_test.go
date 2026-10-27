package fast

import (
	"encoding/hex"
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func fixedSigmaBootstrapNoiseTestParameters(t *testing.T) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39, 40, 39, 40, 60, 60}, LogP: []int{61}, LogDefaultScale: 45,
	})
	require.NoError(t, err)
	return params
}

func fixedSigmaBootstrapNoiseTestCiphertext(params ckks.Parameters, level int) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, level)
	ct.Scale = rlwe.NewScale(math.Exp2(45))
	ct.IsNTT = true
	return ct
}

func fixedSigmaBootstrapNoiseLogN16TestParameters(t *testing.T) ckks.Parameters {
	return fixedSigmaBootstrapNoiseParametersForLogN(t, 16)
}

func fixedSigmaBootstrapNoiseParametersForLogN(t *testing.T, logN int) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: logN,
		Q: []uint64{
			36028797019488257, 549754109953, 549753978881, 549753716737, 549753192449,
			1152921504606584833, 1152921504614055937, 1152921504598720513, 1152921504615628801,
			1152921504616808449, 1152921504597016577, 1152921504595968001, 1152921504618381313,
			72057594038321153, 72057594036879361, 72057594035306497, 72057594040680449,
		},
		P:               []uint64{2305843009211596801, 2305843009210023937, 2305843009218281473, 2305843009208713217, 2305843009218936833},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	require.Equal(t, 1<<logN, params.N())
	require.Len(t, params.Q(), 17)
	require.Len(t, params.P(), 5)
	return params
}

func fixedSigmaBootstrapNoiseLogN16Input(params ckks.Parameters) *rlwe.Ciphertext {
	const level = 12
	ct := ckks.NewCiphertext(params, 1, level)
	ct.Scale = rlwe.NewScale(math.Exp2(60))
	ct.IsNTT, ct.IsMontgomery, ct.IsBatched = true, true, true
	for row := 0; row < 4; row++ {
		for i := range ct.Value[0].Coeffs[row] {
			value := int64(i%3 - 1)
			q := params.RingQ().SubRings[row]
			if value < 0 {
				ct.Value[0].Coeffs[row][i] = q.Modulus - uint64(-value)
			} else {
				ct.Value[0].Coeffs[row][i] = uint64(value)
			}
		}
		params.RingQ().SubRings[row].NTT(ct.Value[0].Coeffs[row], ct.Value[0].Coeffs[row])
		params.RingQ().SubRings[row].MForm(ct.Value[0].Coeffs[row], ct.Value[0].Coeffs[row])
	}
	return ct
}

func TestFixedSigmaBootstrapNoiseLogN16MultiplyIsIndependentAndAliasSafe(t *testing.T) {
	params := fixedSigmaBootstrapNoiseLogN16TestParameters(t)
	initial := fixedSigmaBootstrapNoiseLogN16Input(params)
	expectedScale := rlwe.NewScale(math.Exp2(60))
	rows, err := QPrefixWidth(initial.Level())
	require.NoError(t, err)
	require.Equal(t, 4, rows)
	run := func(aliasOutput bool) (*rlwe.Ciphertext, FixedSigmaBootstrapNoiseSnapshot, *rlwe.Ciphertext) {
		eval := NewEvaluator(params)
		require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 20261027))
		left := initial.CopyNew()
		before := left.CopyNew()
		out := ckks.NewCiphertext(params, 1, left.Level())
		out.IsNTT, out.IsMontgomery, out.IsBatched = true, true, true
		if aliasOutput {
			out = left
		}
		require.NoError(t, eval.MulRelinElementQPrefixRows(left, left.El(), rows, out))
		return out, eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot(), before
	}
	want, wantEvidence, originalBefore := run(false)
	got, gotEvidence, aliasInputBefore := run(true)
	require.Equal(t, want, got, "in-place square must operate on independent noisy operand copies")
	require.Equal(t, originalBefore, initial, "the caller's original input must remain untouched")
	require.Equal(t, aliasInputBefore, initial, "the aliased square must not mutate the original source")
	for _, snapshot := range []FixedSigmaBootstrapNoiseSnapshot{wantEvidence, gotEvidence} {
		require.True(t, snapshot.Enabled)
		require.Equal(t, uint64(1), snapshot.CtCtMultiplyCalls)
		require.Equal(t, uint64(2), snapshot.InjectedPolynomials)
		require.Len(t, snapshot.Events, 1)
		event := snapshot.Events[0]
		require.True(t, event.MultiplyCompleted)
		require.Equal(t, 12, event.Level)
		require.Equal(t, expectedScale.Value.Text('g', -1), event.Scale)
		require.NotEqual(t, event.Left.Seed, event.Right.Seed)
		require.NotEqual(t, event.Left.CoefficientSHA256, event.Right.CoefficientSHA256)
		require.Equal(t, 4, event.Left.Capacity.Rows)
		require.Equal(t, 4, event.Right.Capacity.Rows)
		require.True(t, event.Left.Capacity.StrictFit)
		require.True(t, event.Right.Capacity.StrictFit)
	}
	leftEvent, rightEvent := wantEvidence.Events[0], gotEvidence.Events[0]
	leftEvent.Left.Phases, leftEvent.Right.Phases = nil, nil
	rightEvent.Left.Phases, rightEvent.Right.Phases = nil, nil
	require.Equal(t, leftEvent, rightEvent)
}

func TestFixedSigmaBootstrapNoiseLogN16DefaultOffIsBitExact(t *testing.T) {
	params := fixedSigmaBootstrapNoiseLogN16TestParameters(t)
	input0 := fixedSigmaBootstrapNoiseLogN16Input(params)
	input1 := input0.CopyNew()
	run := func(eval *Evaluator, input *rlwe.Ciphertext) *rlwe.Ciphertext {
		out := ckks.NewCiphertext(params, 1, input.Level())
		out.IsNTT, out.IsMontgomery, out.IsBatched = true, true, true
		require.NoError(t, eval.MulRelinElementQPrefixRows(input, input.El(), 4, out))
		return out
	}
	baseline := NewEvaluator(params)
	defaultOff := NewEvaluator(params)
	require.NoError(t, defaultOff.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 20261027))
	defaultOff.DisableFixedSigmaBootstrapNoiseDiagnostic()
	want := run(baseline, input0)
	got := run(defaultOff, input1)
	require.Equal(t, want, got)
	require.False(t, defaultOff.FixedSigmaBootstrapNoiseDiagnosticSnapshot().Enabled)
}

func TestFixedSigmaBootstrapNoiseLogN16LeavesOtherOperationsUntouched(t *testing.T) {
	params := fixedSigmaBootstrapNoiseLogN16TestParameters(t)
	eval := NewEvaluator(params)
	require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 20261027))
	ct := fixedSigmaBootstrapNoiseLogN16Input(params)
	plain := rlwe.NewElement(eval.GetRLWEParameters(), 0, ct.Level())
	plain.IsNTT, plain.IsMontgomery, plain.IsBatched = true, true, true
	plain.Scale = rlwe.NewScale(1)
	prod := ckks.NewCiphertext(params, 1, ct.Level())
	prod.IsNTT, prod.IsMontgomery, prod.IsBatched = true, true, true
	require.NoError(t, eval.MulElementQPrefixRows(ct, plain, 4, prod))
	scalarOut := ckks.NewCiphertext(params, 1, ct.Level())
	scalarOut.IsNTT, scalarOut.IsMontgomery, scalarOut.IsBatched = true, true, true
	require.NoError(t, eval.Mul(ct, int64(7), scalarOut))
	rotationOut := ckks.NewCiphertext(params, 1, ct.Level())
	rotationOut.IsNTT, rotationOut.IsMontgomery, rotationOut.IsBatched = true, true, true
	require.NoError(t, eval.Automorphism(ct, rotationOut, params.GaloisElementForRotation(1)))
	snapshot := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
	require.Zero(t, snapshot.CtCtMultiplyCalls)
	require.Zero(t, snapshot.InjectedPolynomials)
	require.Empty(t, snapshot.Events)
}

func TestFixedSigmaBootstrapNoiseSupportsLogN13Through16WithEncodedMontgomeryInputs(t *testing.T) {
	const (
		level = 12
		sigma = 4.004280012097432
		seed  = int64(20261027)
	)
	for _, logN := range []int{13, 14, 15, 16} {
		t.Run(fmt.Sprintf("logn%d", logN), func(t *testing.T) {
			params := fixedSigmaBootstrapNoiseParametersForLogN(t, logN)
			encoder := ckks.NewEncoder(params)
			plain := ckks.NewPlaintext(params, level)
			plain.Scale = rlwe.NewScale(math.Exp2(60))
			plain.IsNTT = true
			values := make([]complex128, params.MaxDimensions().Cols)
			values[0], values[1] = complex(0.125, -0.25), complex(-0.0625, 0.03125)
			require.NoError(t, encoder.Encode(values, plain))
			input := ckks.NewCiphertext(params, 1, level)
			input.Value[0] = *plain.Value.CopyNew()
			input.Scale = plain.Scale
			input.IsNTT, input.IsBatched = plain.IsNTT, plain.IsBatched
			input.LogDimensions, input.IsBitReversed = plain.LogDimensions, plain.IsBitReversed
			input.IsMontgomery = true
			for row := 0; row < 4; row++ {
				params.RingQ().SubRings[row].MForm(input.Value[0].Coeffs[row], input.Value[0].Coeffs[row])
			}
			original := input.CopyNew()

			// A real encoded zero-secret c0 traverses ct×pt, scalar and rotation
			// operations without triggering the ct×ct-only hook.
			untouchedEval := NewEvaluator(params)
			require.NoError(t, untouchedEval.EnableFixedSigmaBootstrapNoiseDiagnostic(sigma, seed))
			plainOperand := rlwe.NewElement(untouchedEval.GetRLWEParameters(), 0, level)
			plainOperand.IsNTT, plainOperand.IsMontgomery, plainOperand.IsBatched = true, true, true
			plainOperand.Scale = rlwe.NewScale(1)
			ctPtOut := ckks.NewCiphertext(params, 1, level)
			ctPtOut.IsNTT, ctPtOut.IsMontgomery, ctPtOut.IsBatched = true, true, true
			require.NoError(t, untouchedEval.MulElementQPrefixRows(input, plainOperand, 4, ctPtOut))
			scalarOut := ckks.NewCiphertext(params, 1, level)
			scalarOut.IsNTT, scalarOut.IsMontgomery, scalarOut.IsBatched = true, true, true
			require.NoError(t, untouchedEval.Mul(input, int64(7), scalarOut))
			rotationOut := ckks.NewCiphertext(params, 1, level)
			rotationOut.IsNTT, rotationOut.IsMontgomery, rotationOut.IsBatched = true, true, true
			require.NoError(t, untouchedEval.Automorphism(input, rotationOut, params.GaloisElementForRotation(1)))
			snapshot := untouchedEval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
			require.Zero(t, snapshot.CtCtMultiplyCalls)
			require.Zero(t, snapshot.InjectedPolynomials)
			require.Empty(t, snapshot.Events)

			// Disabling the opt-in hook must preserve the ordinary Fast product
			// bit-for-bit at every tested ring dimension.
			baseline := NewEvaluator(params)
			disabled := NewEvaluator(params)
			require.NoError(t, disabled.EnableFixedSigmaBootstrapNoiseDiagnostic(sigma, seed))
			disabled.DisableFixedSigmaBootstrapNoiseDiagnostic()
			baselineOut := ckks.NewCiphertext(params, 1, level)
			disabledOut := ckks.NewCiphertext(params, 1, level)
			for _, out := range []*rlwe.Ciphertext{baselineOut, disabledOut} {
				out.IsNTT, out.IsMontgomery, out.IsBatched = true, true, true
			}
			rows, err := QPrefixWidth(level)
			require.NoError(t, err)
			require.NoError(t, baseline.MulRelinElementQPrefixRows(input.CopyNew(), input.El(), rows, baselineOut))
			require.NoError(t, disabled.MulRelinElementQPrefixRows(input.CopyNew(), input.El(), rows, disabledOut))
			require.Equal(t, baselineOut, disabledOut)
			require.False(t, disabled.FixedSigmaBootstrapNoiseDiagnosticSnapshot().Enabled)

			runSquare := func(aliasOutput bool) (*rlwe.Ciphertext, FixedSigmaBootstrapNoiseSnapshot) {
				eval := NewEvaluator(params)
				require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(sigma, seed))
				left := input.CopyNew()
				out := ckks.NewCiphertext(params, 1, level)
				out.IsNTT, out.IsMontgomery, out.IsBatched = true, true, true
				if aliasOutput {
					out = left
				}
				require.NoError(t, eval.MulRelinElementQPrefixRows(left, left.El(), rows, out))
				return out, eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
			}
			copiedOut, copiedEvidence := runSquare(false)
			aliasedOut, aliasedEvidence := runSquare(true)
			require.Equal(t, copiedOut, aliasedOut, "square output alias must use independent operand copies")
			require.Equal(t, original, input, "hook must not mutate original encoded input")
			for _, result := range []struct {
				name     string
				observed FixedSigmaBootstrapNoiseSnapshot
			}{{"copied", copiedEvidence}, {"aliased", aliasedEvidence}} {
				require.True(t, result.observed.Enabled, result.name)
				require.Equal(t, uint64(1), result.observed.CtCtMultiplyCalls, result.name)
				require.Equal(t, uint64(2), result.observed.InjectedPolynomials, result.name)
				require.Len(t, result.observed.Events, 1, result.name)
				event := result.observed.Events[0]
				require.True(t, event.MultiplyCompleted, result.name)
				require.Equal(t, params.N(), event.Left.SampleCount, result.name)
				require.Equal(t, params.N(), event.Right.SampleCount, result.name)
				require.Equal(t, level, event.Level, result.name)
				require.Equal(t, input.Scale.Value.Text('g', -1), event.Scale, result.name)
				require.NotEqual(t, event.Left.Seed, event.Right.Seed, result.name)
				require.NotEqual(t, event.Left.CoefficientSHA256, event.Right.CoefficientSHA256, result.name)
				require.Len(t, mustDecodeFixedSigmaSHA256(t, event.Left.CoefficientSHA256), 32, result.name)
				require.Len(t, mustDecodeFixedSigmaSHA256(t, event.Right.CoefficientSHA256), 32, result.name)
				require.True(t, event.Left.Capacity.StrictFit, result.name)
				require.True(t, event.Right.Capacity.StrictFit, result.name)
				require.Equal(t, rows, event.Left.Capacity.Rows, result.name)
				require.Equal(t, rows, event.Right.Capacity.Rows, result.name)
			}
			require.Equal(t, copiedEvidence.Events[0].Left.CoefficientSHA256, aliasedEvidence.Events[0].Left.CoefficientSHA256)
			require.Equal(t, copiedEvidence.Events[0].Right.CoefficientSHA256, aliasedEvidence.Events[0].Right.CoefficientSHA256)
		})
	}
}

func mustDecodeFixedSigmaSHA256(t *testing.T, value string) []byte {
	t.Helper()
	digest, err := hex.DecodeString(value)
	require.NoError(t, err)
	return digest
}

func TestFixedSigmaBootstrapNoiseRejectsUninitializedParameters(t *testing.T) {
	eval := &Evaluator{}
	err := eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 1)
	require.ErrorContains(t, err, "do not contain an initialized CKKS ring")
	require.False(t, eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot().Enabled)
}

func TestFixedSigmaBootstrapNoiseDefaultOffIsBitExact(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	left, right := fixedSigmaBootstrapNoiseTestCiphertext(params, 3), fixedSigmaBootstrapNoiseTestCiphertext(params, 3)
	rows, err := QPrefixWidth(3)
	require.NoError(t, err)
	want, got := ckks.NewCiphertext(params, 1, 3), ckks.NewCiphertext(params, 1, 3)
	want.IsNTT, got.IsNTT = true, true
	baseline := NewEvaluator(params)
	actual := NewEvaluator(params)
	require.NoError(t, baseline.MulRelinElementQPrefixRows(left, right.El(), rows, want))
	require.NoError(t, actual.MulRelinElementQPrefixRows(left, right.El(), rows, got))
	require.Equal(t, want, got)
	require.False(t, actual.FixedSigmaBootstrapNoiseDiagnosticSnapshot().Enabled)
}

func TestFixedSigmaBootstrapNoiseInjectsIndependentSidesAndPreservesAliasedInputs(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	eval := NewEvaluator(params)
	require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 20261027))
	left := fixedSigmaBootstrapNoiseTestCiphertext(params, 3)
	before := left.CopyNew()
	out := ckks.NewCiphertext(params, 1, 3)
	out.IsNTT = true
	rows, err := QPrefixWidth(3)
	require.NoError(t, err)
	require.NoError(t, eval.MulRelinElementQPrefixRows(left, left.El(), rows, out))
	require.Equal(t, before, left, "the square input must remain unchanged when output does not alias it")

	snapshot := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
	require.True(t, snapshot.Enabled)
	require.Equal(t, 4.004280012097432, snapshot.Sigma)
	require.Equal(t, uint64(1), snapshot.CtCtMultiplyCalls)
	require.Equal(t, uint64(2), snapshot.InjectedPolynomials)
	require.Len(t, snapshot.Events, 1)
	event := snapshot.Events[0]
	require.True(t, event.MultiplyCompleted)
	require.NotEqual(t, event.Left.Seed, event.Right.Seed)
	require.NotEqual(t, event.Left.CoefficientSHA256, event.Right.CoefficientSHA256)
	require.Equal(t, event.Left.Seed, ckks.FixedSigmaBootstrapNoiseSeed(20261027, 1, "left"))
	require.Equal(t, event.Right.Seed, ckks.FixedSigmaBootstrapNoiseSeed(20261027, 1, "right"))
	require.Equal(t, uint64(1), event.Operation)
	require.Equal(t, 4.004280012097432, event.Left.Sigma)
	require.Equal(t, 4.004280012097432, event.Right.Sigma)
	for side, evidence := range map[string]ckks.FixedSigmaCoefficientNoiseEvidence{"left": event.Left, "right": event.Right} {
		for _, phase := range evidence.Phases {
			require.False(t, phase.MemoryStatsCollected, "%s Bootstrap injection must omit per-phase MemStats", side)
			require.Zero(t, phase.AllocBytes)
			require.Zero(t, phase.Allocs)
		}
	}
}

func TestFixedSigmaBootstrapNoiseAliasedSquareOutputMatchesCopiedSquare(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	initial := fixedSigmaBootstrapNoiseTestCiphertext(params, 3)
	rows, err := QPrefixWidth(3)
	require.NoError(t, err)
	run := func(aliasOutput bool) (*rlwe.Ciphertext, FixedSigmaBootstrapNoiseSnapshot) {
		eval := NewEvaluator(params)
		require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 8181))
		left := initial.CopyNew()
		out := ckks.NewCiphertext(params, 1, 3)
		out.IsNTT = true
		if aliasOutput {
			out = left
		}
		require.NoError(t, eval.MulRelinElementQPrefixRows(left, left.El(), rows, out))
		return out, eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
	}
	copiedOutput, copiedEvidence := run(false)
	aliasedOutput, aliasedEvidence := run(true)
	require.Equal(t, copiedOutput, aliasedOutput, "in-place square must use independent noisy operand copies before writing its result")
	require.Len(t, copiedEvidence.Events, 1)
	require.Len(t, aliasedEvidence.Events, 1)
	copiedEvent, aliasedEvent := copiedEvidence.Events[0], aliasedEvidence.Events[0]
	copiedEvent.Left.Phases, copiedEvent.Right.Phases = nil, nil
	aliasedEvent.Left.Phases, aliasedEvent.Right.Phases = nil, nil
	require.Equal(t, copiedEvent, aliasedEvent)
	require.Equal(t, uint64(1), aliasedEvidence.CtCtMultiplyCalls)
}

func TestFixedSigmaBootstrapNoiseLeavesCiphertextPlaintextProductUntouched(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	eval := NewEvaluator(params)
	require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 7))
	ct := fixedSigmaBootstrapNoiseTestCiphertext(params, 3)
	plain := rlwe.NewElement(eval.GetRLWEParameters(), 0, 3)
	plain.IsNTT = true
	plain.IsBatched = ct.IsBatched
	plain.Scale = rlwe.NewScale(1)
	out := ckks.NewCiphertext(params, 1, 3)
	out.IsNTT = true
	require.NoError(t, eval.MulRelinElementQPrefixRows(ct, plain, 4, out))
	snapshot := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
	require.Zero(t, snapshot.CtCtMultiplyCalls)
	require.Zero(t, snapshot.InjectedPolynomials)
	require.Empty(t, snapshot.Events)
}

func TestFixedSigmaBootstrapNoiseHooksLazyMultiplyAndLeavesScalarAlone(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	eval := NewEvaluator(params)
	require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 107))
	left, right := fixedSigmaBootstrapNoiseTestCiphertext(params, 3), fixedSigmaBootstrapNoiseTestCiphertext(params, 3)
	product := ckks.NewCiphertext(params, 2, 3)
	product.IsNTT = true
	require.NoError(t, eval.MulElementQPrefixRows(left, right.El(), 4, product))
	require.Equal(t, 2, product.Degree())
	require.Equal(t, uint64(1), eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot().CtCtMultiplyCalls)

	scalarOut := ckks.NewCiphertext(params, 1, 3)
	scalarOut.IsNTT = true
	require.NoError(t, eval.Mul(left, int64(7), scalarOut))
	snapshot := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
	require.Equal(t, uint64(1), snapshot.CtCtMultiplyCalls)
	require.Equal(t, uint64(2), snapshot.InjectedPolynomials)
	require.Len(t, snapshot.Events, 1)
}

func TestFixedSigmaBootstrapNoiseLeavesRotationUntouched(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	eval := NewEvaluator(params)
	require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 108))
	in := fixedSigmaBootstrapNoiseTestCiphertext(params, 3)
	out := ckks.NewCiphertext(params, 1, 3)
	out.IsNTT = true
	require.NoError(t, eval.Automorphism(in, out, params.GaloisElementForRotation(1)))
	snapshot := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
	require.Zero(t, snapshot.CtCtMultiplyCalls)
	require.Zero(t, snapshot.InjectedPolynomials)
	require.Empty(t, snapshot.Events)
}

func TestFixedSigmaBootstrapNoiseReproducibleAcrossEvaluatorsAndOffToggle(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	left, right := fixedSigmaBootstrapNoiseTestCiphertext(params, 2), fixedSigmaBootstrapNoiseTestCiphertext(params, 2)
	run := func(enabled bool) (*rlwe.Ciphertext, FixedSigmaBootstrapNoiseSnapshot) {
		eval := NewEvaluator(params)
		if enabled {
			require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, -91))
		}
		out := ckks.NewCiphertext(params, 1, 2)
		out.IsNTT = true
		require.NoError(t, eval.MulRelinElementQPrefixRows(left.CopyNew(), right.CopyNew().El(), 3, out))
		return out, eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
	}
	first, firstEvidence := run(true)
	second, secondEvidence := run(true)
	require.Equal(t, first, second)
	require.Equal(t, firstEvidence.Enabled, secondEvidence.Enabled)
	require.Equal(t, firstEvidence.Sigma, secondEvidence.Sigma)
	require.Equal(t, firstEvidence.RootSeed, secondEvidence.RootSeed)
	require.Equal(t, firstEvidence.CtCtMultiplyCalls, secondEvidence.CtCtMultiplyCalls)
	require.Equal(t, firstEvidence.InjectedPolynomials, secondEvidence.InjectedPolynomials)
	require.Len(t, firstEvidence.Events, 1)
	require.Len(t, secondEvidence.Events, 1)
	firstEvent, secondEvent := firstEvidence.Events[0], secondEvidence.Events[0]
	firstEvent.Left.Phases, firstEvent.Right.Phases = nil, nil
	secondEvent.Left.Phases, secondEvent.Right.Phases = nil, nil
	require.Equal(t, firstEvent, secondEvent)
	require.Equal(t, 1, len(firstEvidence.Events))
	base, disabled := run(false)
	require.False(t, disabled.Enabled)
	require.NotEqual(t, first, base)

	zeroEvaluator := NewEvaluator(params)
	require.Error(t, zeroEvaluator.EnableFixedSigmaBootstrapNoiseDiagnostic(0, 1))
	require.False(t, zeroEvaluator.FixedSigmaBootstrapNoiseDiagnosticSnapshot().Enabled)
}

func TestFixedSigmaBootstrapNoiseUsesActiveRowsAcrossLevels(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	for _, level := range []int{1, 2, 3, 4} {
		t.Run(string(rune('0'+level)), func(t *testing.T) {
			eval := NewEvaluator(params)
			require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 55))
			left, right := fixedSigmaBootstrapNoiseTestCiphertext(params, level), fixedSigmaBootstrapNoiseTestCiphertext(params, level)
			out := ckks.NewCiphertext(params, 1, level)
			out.IsNTT = true
			rows, err := QPrefixWidth(level)
			require.NoError(t, err)
			require.NoError(t, eval.MulRelinElementQPrefixRows(left, right.El(), rows, out))
			evidence := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
			require.Equal(t, rows, evidence.Events[0].Left.Capacity.Rows)
			require.Equal(t, rows, evidence.Events[0].Right.Capacity.Rows)
		})
	}
}

func TestFixedSigmaBootstrapNoiseSupportsMontgomeryWithoutMutatingOperands(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	eval := NewEvaluator(params)
	require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 4))
	left, right := fixedSigmaBootstrapNoiseTestCiphertext(params, 3), fixedSigmaBootstrapNoiseTestCiphertext(params, 3)
	for _, ct := range []*rlwe.Ciphertext{left, right} {
		for d := range ct.Value {
			for row := 0; row < 4; row++ {
				params.RingQ().SubRings[row].MForm(ct.Value[d].Coeffs[row], ct.Value[d].Coeffs[row])
			}
		}
		ct.IsMontgomery = true
	}
	leftBefore, rightBefore := left.CopyNew(), right.CopyNew()
	out := ckks.NewCiphertext(params, 1, 3)
	out.IsNTT, out.IsMontgomery = true, true
	require.NoError(t, eval.MulRelinElementQPrefixRows(left, right.El(), 4, out))
	require.Equal(t, leftBefore, left)
	require.Equal(t, rightBefore, right)
	snapshot := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot()
	require.Equal(t, uint64(1), snapshot.CtCtMultiplyCalls)
	require.Equal(t, uint64(2), snapshot.InjectedPolynomials)
	require.Len(t, snapshot.Events, 1)
	require.True(t, snapshot.Events[0].MultiplyCompleted)
	require.Contains(t, snapshot.Events[0].Left.DispatchEvidence, "MForm-mapped")
	require.Contains(t, snapshot.Events[0].Right.DispatchEvidence, "MForm-mapped")
}
