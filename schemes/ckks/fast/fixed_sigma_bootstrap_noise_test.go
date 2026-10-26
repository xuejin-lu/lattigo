package fast

import (
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
	require.Equal(t, event.Left.Seed, fixedSigmaBootstrapNoiseSeed(20261027, 1, "left"))
	require.Equal(t, event.Right.Seed, fixedSigmaBootstrapNoiseSeed(20261027, 1, "right"))
	require.Equal(t, uint64(1), event.Operation)
	require.Equal(t, 4.004280012097432, event.Left.Sigma)
	require.Equal(t, 4.004280012097432, event.Right.Sigma)
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

func TestFixedSigmaBootstrapNoiseRejectsMontgomeryWithoutMutatingOperands(t *testing.T) {
	params := fixedSigmaBootstrapNoiseTestParameters(t)
	eval := NewEvaluator(params)
	require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 4))
	left, right := fixedSigmaBootstrapNoiseTestCiphertext(params, 3), fixedSigmaBootstrapNoiseTestCiphertext(params, 3)
	for d := range left.Value {
		for row := 0; row < 4; row++ {
			params.RingQ().SubRings[row].MForm(left.Value[d].Coeffs[row], left.Value[d].Coeffs[row])
		}
	}
	left.IsMontgomery = true
	leftBefore, rightBefore := left.CopyNew(), right.CopyNew()
	out := ckks.NewCiphertext(params, 1, 3)
	out.IsNTT = true
	require.ErrorContains(t, eval.MulRelinElementQPrefixRows(left, right.El(), 4, out), "Montgomery")
	require.Equal(t, leftBefore, left)
	require.Equal(t, rightBefore, right)
	require.Zero(t, eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot().CtCtMultiplyCalls)
}
