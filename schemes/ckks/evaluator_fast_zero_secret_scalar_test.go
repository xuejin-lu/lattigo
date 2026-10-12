package ckks

import (
	"crypto/sha256"
	"encoding/binary"
	"math"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

type scalarCoreCall struct {
	name string
	rows int
}

type recordingScalarCore struct {
	delegate fastcore.ScalarCore
	calls    []scalarCoreCall
}

func (recorder *recordingScalarCore) AddSubRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows int, scalar *bignum.Complex, sub bool) error {
	recorder.calls = append(recorder.calls, scalarCoreCall{name: "addsub", rows: rows})
	return recorder.delegate.AddSubRows(ringQ, op0, opOut, rows, scalar, sub)
}

func (recorder *recordingScalarCore) MulRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows, levelsConsumed int, scalar *bignum.Complex) error {
	recorder.calls = append(recorder.calls, scalarCoreCall{name: "mul", rows: rows})
	return recorder.delegate.MulRows(ringQ, op0, opOut, rows, levelsConsumed, scalar)
}

func (recorder *recordingScalarCore) MulAtScaleRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows int, scalar *bignum.Complex, scale rlwe.Scale) error {
	recorder.calls = append(recorder.calls, scalarCoreCall{name: "mul_at_scale", rows: rows})
	return recorder.delegate.MulAtScaleRows(ringQ, op0, opOut, rows, scalar, scale)
}

func (recorder *recordingScalarCore) MulThenAddRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows, levelsConsumed int, scalar *bignum.Complex) error {
	recorder.calls = append(recorder.calls, scalarCoreCall{name: "mul_then_add", rows: rows})
	return recorder.delegate.MulThenAddRows(ringQ, op0, opOut, rows, levelsConsumed, scalar)
}

func (recorder *recordingScalarCore) MulThenAddEncodedRows(ringQ *ring.Ring, source []ring.Poly, opOut *rlwe.Ciphertext, rows int, values [fastcore.MaxQPrefixWidth][2]uint64) error {
	recorder.calls = append(recorder.calls, scalarCoreCall{name: "mul_then_add_encoded", rows: rows})
	return recorder.delegate.MulThenAddEncodedRows(ringQ, source, opOut, rows, values)
}

func btsCoverageLogN16Parameters(t *testing.T) Parameters {
	t.Helper()
	logQ := []int{42}
	for i := 0; i < 3; i++ {
		logQ = append(logQ, 32)
	}
	for i := 0; i < 15; i++ {
		logQ = append(logQ, 32)
	}
	for i := 0; i < 8; i++ {
		logQ = append(logQ, 42)
	}
	for i := 0; i < 3; i++ {
		logQ = append(logQ, 42)
	}
	logP := []int{42, 42, 42, 42, 42}
	params, err := NewParametersFromLiteral(ParametersLiteral{LogN: 16, LogQ: logQ, LogP: logP, LogDefaultScale: 32})
	require.NoError(t, err)
	require.Len(t, params.Q(), 30)
	require.Equal(t, 29, params.MaxLevel())
	return params
}

func btsCoverageCiphertext(t *testing.T, params Parameters, level int, seed uint64, fullBacking bool) *rlwe.Ciphertext {
	t.Helper()
	ct := fastcore.NewCompactCiphertext(params, 1, level)
	if fullBacking {
		ct = NewCiphertext(params, 1, level)
	}
	ct.IsNTT = true
	ct.IsMontgomery = false
	ct.IsBatched = true
	ct.Scale = rlwe.NewScale(math.Exp2(32))
	ct.LogDimensions = ring.Dimensions{Cols: params.LogMaxSlots()}
	rows, err := fastcore.QPrefixWidth(level)
	require.NoError(t, err)
	for row := 0; row < rows; row++ {
		modulus := params.RingQ().SubRings[row].Modulus
		for i := range ct.Value[0].Coeffs[row] {
			ct.Value[0].Coeffs[row][i] = (seed + uint64(19*row+7*i+1)) % modulus
		}
		clear(ct.Value[1].Coeffs[row])
	}
	if fullBacking {
		for component := range ct.Value {
			for row := rows; row <= level; row++ {
				for i := range ct.Value[component].Coeffs[row] {
					ct.Value[component].Coeffs[row][i] = ^uint64(0) - seed - uint64(31*component+13*row+i)
				}
			}
		}
	}
	return ct
}

func btsCoverageScalarOracle(t *testing.T, params Parameters, ct *rlwe.Ciphertext, scalar rlwe.Operand, scale rlwe.Scale, op string) []ring.Poly {
	t.Helper()
	rows, err := fastcore.QPrefixWidth(ct.Level())
	require.NoError(t, err)
	ringQ := params.RingQ().AtLevel(rows - 1)
	rnsReal, rnsImag := bigComplexToRNSScalar(ringQ, &scale.Value, bignum.ToComplex(scalar, params.EncodingPrecision()))
	// This is the Standard evaluator's complex-scalar embedding step, followed
	// by the independent generic ring primitive used as the oracle.
	for i, subring := range ringQ.SubRings[:rows] {
		rnsImag[i] = ring.MRed(rnsImag[i], subring.RootsForward[1], subring.Modulus, subring.MRedConstant)
		rnsReal[i], rnsImag[i] = ring.CRed(rnsReal[i]+rnsImag[i], subring.Modulus), ring.CRed(rnsReal[i]+subring.Modulus-rnsImag[i], subring.Modulus)
	}
	want := make([]ring.Poly, len(ct.Value))
	for component := range ct.Value {
		source := ring.NewPoly(params.N(), rows-1)
		want[component] = ring.NewPoly(params.N(), rows-1)
		for row := 0; row < rows; row++ {
			copy(source.Coeffs[row], ct.Value[component].Coeffs[row])
		}
		switch op {
		case "mul":
			ringQ.MulDoubleRNSScalar(source, rnsReal, rnsImag, want[component])
		case "add":
			if component == 0 {
				ringQ.AddDoubleRNSScalar(source, rnsReal, rnsImag, want[component])
			} else {
				want[component].Copy(source)
			}
		case "sub":
			if component == 0 {
				ringQ.SubDoubleRNSScalar(source, rnsReal, rnsImag, want[component])
			} else {
				want[component].Copy(source)
			}
		default:
			t.Fatalf("unknown scalar oracle operation %q", op)
		}
	}
	return want
}

// fastCKKSQPrefixDecodeView makes the authoritative compact rows decryptable
// at their own logical level for generic CKKS accuracy tests. It never fills
// or reads dormant rows.
func fastCKKSQPrefixDecodeView(t *testing.T, params Parameters, ct *rlwe.Ciphertext) *rlwe.Ciphertext {
	t.Helper()
	rows, err := fastcore.QPrefixWidth(ct.Level())
	require.NoError(t, err)
	view := NewCiphertext(params, ct.Degree(), rows-1)
	*view.MetaData = *ct.MetaData
	for component := range ct.Value {
		for row := 0; row < rows; row++ {
			copy(view.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
		}
	}
	return view
}

func requireBTSCoverageRows(t *testing.T, params Parameters, got *rlwe.Ciphertext, want []ring.Poly) {
	t.Helper()
	rows, err := fastcore.QPrefixWidth(got.Level())
	require.NoError(t, err)
	requireBTSActiveCoverageRows(t, got, rows, want)
	for component := range got.Value {
		for row := rows; row <= got.Level(); row++ {
			require.Nil(t, got.Value[component].Coeffs[row], "dormant component=%d q%d", component, row)
		}
	}
}

func requireBTSActiveCoverageRows(t *testing.T, got *rlwe.Ciphertext, rows int, want []ring.Poly) {
	t.Helper()
	require.Len(t, got.Value, len(want))
	for component := range want {
		for row := 0; row < rows; row++ {
			require.Len(t, got.Value[component].Coeffs[row], len(want[component].Coeffs[row]), "component=%d q%d", component, row)
			for i, expected := range want[component].Coeffs[row] {
				require.Equal(t, expected, got.Value[component].Coeffs[row][i], "component=%d q%d coefficient=%d", component, row, i)
			}
		}
	}
}

func btsDormantRowsDigest(ct *rlwe.Ciphertext, rows int) [32]byte {
	h := sha256.New()
	var word [8]byte
	for component := range ct.Value {
		for row := rows; row < len(ct.Value[component].Coeffs); row++ {
			coeffs := ct.Value[component].Coeffs[row]
			if coeffs == nil {
				_, _ = h.Write([]byte{0})
				continue
			}
			_, _ = h.Write([]byte{1})
			for _, coefficient := range coeffs {
				binary.LittleEndian.PutUint64(word[:], coefficient)
				_, _ = h.Write(word[:])
			}
		}
	}
	var digest [32]byte
	copy(digest[:], h.Sum(nil))
	return digest
}

func requireBTSCoverageCiphertextPayloadEqual(t *testing.T, expected, got *rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, len(expected.Value), len(got.Value))
	require.Equal(t, expected.Level(), got.Level())
	require.Equal(t, expected.IsNTT, got.IsNTT)
	require.Equal(t, expected.IsMontgomery, got.IsMontgomery)
	require.Equal(t, expected.IsBatched, got.IsBatched)
	require.Equal(t, expected.MetaData, got.MetaData)
	require.True(t, expected.Scale.Equal(got.Scale))
	for component := range expected.Value {
		require.Equal(t, len(expected.Value[component].Coeffs), len(got.Value[component].Coeffs))
		for row := range expected.Value[component].Coeffs {
			wantRow, gotRow := expected.Value[component].Coeffs[row], got.Value[component].Coeffs[row]
			if wantRow == nil {
				require.Nil(t, gotRow, "component=%d q%d backing", component, row)
				continue
			}
			require.Len(t, gotRow, len(wantRow), "component=%d q%d", component, row)
			for i, value := range wantRow {
				require.Equal(t, value, gotRow[i], "component=%d q%d coefficient=%d", component, row, i)
			}
		}
	}
}

func TestPublicFastBTSLogN16ScalarAPIsUseSharedQPrefixAndPreserveCallChain(t *testing.T) {
	params := btsCoverageLogN16Parameters(t)
	level := params.MaxLevel() - 3
	rows, err := fastcore.QPrefixWidth(level)
	require.NoError(t, err)
	require.Equal(t, 4, rows)
	evaluator := NewEvaluator(params, nil)
	recorder := &recordingScalarCore{delegate: fastcore.NewScalarWorkspace(params.RingQ(), rows)}
	evaluator.scalarCore = recorder

	// The full backing contains deliberately poisoned dormant rows. Mul(-1i)
	// must match the generic ring oracle on q0..q3 and compact the result.
	full := btsCoverageCiphertext(t, params, level, 23, true)
	mulScale := rlwe.NewScale(1) // -1i is a Gaussian integer and consumes no scale.
	wantMinusI := btsCoverageScalarOracle(t, params, full, complex(0, -1), mulScale, "mul")
	originalScale := full.Scale
	require.NoError(t, evaluator.Mul(full, complex(0, -1), full))
	requireBTSCoverageRows(t, params, full, wantMinusI)
	require.Equal(t, level, full.Level())
	require.True(t, originalScale.Equal(full.Scale))
	require.True(t, full.IsNTT)
	require.False(t, full.IsMontgomery)
	for _, coefficient := range full.Value[1].Coeffs[:rows] {
		for _, value := range coefficient {
			require.Zero(t, value)
		}
	}
	require.Equal(t, scalarCoreCall{name: "mul", rows: rows}, recorder.calls[0])

	// A downstream ordinary Add must accept the compact scalar result. BTS
	// computes its real-branch Add separately; the real/imaginary branches are
	// parallel outputs, not sequential checkpoints.
	zV := btsCoverageCiphertext(t, params, level, 97, false)
	zV.Scale = full.Scale
	wantAdd := make([]ring.Poly, len(full.Value))
	for component := range full.Value {
		wantAdd[component] = ring.NewPoly(params.N(), rows-1)
		for row := 0; row < rows; row++ {
			params.RingQ().SubRings[row].Add(full.Value[component].Coeffs[row], zV.Value[component].Coeffs[row], wantAdd[component].Coeffs[row])
		}
	}
	require.NoError(t, evaluator.Add(full, zV, full))
	requireBTSCoverageRows(t, params, full, wantAdd)

	// Add/Sub scalar and MulThenAdd are the other scalar dispatches in EvalMod.
	before := btsCoverageCiphertext(t, params, level, 113, false)
	addWant := btsCoverageScalarOracle(t, params, before, int64(1), before.Scale, "add")
	require.NoError(t, evaluator.Add(before, int64(1), before))
	requireBTSCoverageRows(t, params, before, addWant)
	subWant := btsCoverageScalarOracle(t, params, before, int64(1), before.Scale, "sub")
	require.NoError(t, evaluator.Sub(before, int64(1), before))
	requireBTSCoverageRows(t, params, before, subWant)

	accumulator := btsCoverageCiphertext(t, params, level, 149, false)
	source := btsCoverageCiphertext(t, params, level, 173, false)
	wantProduct := btsCoverageScalarOracle(t, params, source, int64(2), rlwe.NewScale(1), "mul")
	for component := range wantProduct {
		for row := 0; row < rows; row++ {
			params.RingQ().SubRings[row].Add(wantProduct[component].Coeffs[row], accumulator.Value[component].Coeffs[row], wantProduct[component].Coeffs[row])
		}
	}
	accScale := accumulator.Scale
	require.NoError(t, evaluator.MulThenAdd(source, int64(2), accumulator))
	requireBTSCoverageRows(t, params, accumulator, wantProduct)
	require.True(t, accScale.Equal(accumulator.Scale))
	require.Equal(t, scalarCoreCall{name: "mul_then_add", rows: rows}, recorder.calls[3])

	levelBeforeRescale := accumulator.Level()
	require.NoError(t, evaluator.Rescale(accumulator, accumulator))
	require.Equal(t, levelBeforeRescale-params.LevelsConsumedPerRescaling(), accumulator.Level())
	require.Equal(t, []scalarCoreCall{{name: "mul", rows: rows}, {name: "addsub", rows: rows}, {name: "addsub", rows: rows}, {name: "mul_then_add", rows: rows}}, recorder.calls)
}

func TestPublicFastBTSLogN16ComplexScalarMulPreservesMontgomeryRepresentation(t *testing.T) {
	params := btsCoverageLogN16Parameters(t)
	level := params.MaxLevel() - 3
	input := btsCoverageCiphertext(t, params, level, 191, true)
	input.IsMontgomery = true
	rows, err := fastcore.QPrefixWidth(level)
	require.NoError(t, err)
	for component := range input.Value {
		for row := 0; row < rows; row++ {
			subring := params.RingQ().SubRings[row]
			for i, value := range input.Value[component].Coeffs[row] {
				input.Value[component].Coeffs[row][i] = ring.MForm(value, subring.Modulus, subring.BRedConstant)
			}
		}
	}
	mulScale := rlwe.NewScale(1)
	want := btsCoverageScalarOracle(t, params, input, complex(0, -1), mulScale, "mul")
	output := fastcore.NewCompactCiphertext(params, 1, level)
	output.IsNTT, output.IsMontgomery = true, true
	require.NoError(t, NewEvaluator(params, nil).Mul(input, complex(0, -1), output))
	requireBTSCoverageRows(t, params, output, want)
	require.True(t, output.IsNTT)
	require.True(t, output.IsMontgomery)
	require.Equal(t, level, output.Level())
	require.True(t, input.Scale.Equal(output.Scale))
}

func TestPublicFastBTSLogN16MulThenAddInPlaceScalePromotion(t *testing.T) {
	params := btsCoverageLogN16Parameters(t)
	level := params.MaxLevel() - 3
	rows, err := fastcore.QPrefixWidth(level)
	require.NoError(t, err)
	mulScale, err := fastcore.CoefficientScale(params.RingQ(), level, params.LevelsConsumedPerRescaling())
	require.NoError(t, err)
	input := btsCoverageCiphertext(t, params, level, 229, true)
	dormantBefore := btsDormantRowsDigest(input, rows)
	original := btsCoverageScalarOracle(t, params, input, int64(1), rlwe.NewScale(1), "mul")
	term := btsCoverageScalarOracle(t, params, input, complex(0.5, 0), mulScale, "mul")
	ringQ := params.RingQ().AtLevel(rows - 1)
	want := make([]ring.Poly, len(input.Value))
	for component := range input.Value {
		scaledAccumulator := ring.NewPoly(params.N(), rows-1)
		want[component] = ring.NewPoly(params.N(), rows-1)
		ringQ.MulScalarBigint(original[component], mulScale.BigInt(), scaledAccumulator)
		ringQ.Add(scaledAccumulator, term[component], want[component])
	}
	originalScale := input.Scale
	require.NoError(t, NewEvaluator(params, nil).MulThenAdd(input, complex(0.5, 0), input))
	requireBTSActiveCoverageRows(t, input, rows, want)
	require.Equal(t, dormantBefore, btsDormantRowsDigest(input, rows), "MulThenAdd must not read or rewrite dormant rows")
	require.True(t, originalScale.Mul(mulScale).Equal(input.Scale))
	require.Equal(t, level, input.Level())
	require.True(t, input.IsNTT)
	require.False(t, input.IsMontgomery)
	require.NoError(t, NewEvaluator(params, nil).Rescale(input, input))
	require.Equal(t, level-params.LevelsConsumedPerRescaling(), input.Level())
}

func TestPublicFastBTSLogN16ScalarRejectsMissingAuthoritativeRowTransactionally(t *testing.T) {
	params := btsCoverageLogN16Parameters(t)
	level := params.MaxLevel() - 3
	evaluator := NewEvaluator(params, nil)
	input := btsCoverageCiphertext(t, params, level, 211, false)
	output := btsCoverageCiphertext(t, params, level, 239, false)
	wantInput, wantOutput := input.CopyNew(), output.CopyNew()
	for component := range input.Value {
		for row := range input.Value[component].Coeffs {
			if input.Value[component].Coeffs[row] == nil {
				wantInput.Value[component].Coeffs[row] = nil
			}
		}
	}
	for component := range output.Value {
		for row := range output.Value[component].Coeffs {
			if output.Value[component].Coeffs[row] == nil {
				wantOutput.Value[component].Coeffs[row] = nil
			}
		}
	}
	input.Value[0].Coeffs[2] = nil
	wantInput.Value[0].Coeffs[2] = nil
	err := evaluator.Mul(input, complex(0, -1), output)
	require.ErrorContains(t, err, "q2")
	requireBTSCoverageCiphertextPayloadEqual(t, wantInput, input)
	requireBTSCoverageCiphertextPayloadEqual(t, wantOutput, output)

	// Non-scalar MulThenAdd remains a complete-Q compatibility lane. A compact
	// input must fail before its full-Q ring kernel can read dormant rows.
	plaintext := NewPlaintext(params, level)
	fullOutput := NewCiphertext(params, 1, level)
	fullOutputBefore := fullOutput.CopyNew()
	err = evaluator.MulThenAdd(input, plaintext, fullOutput)
	require.ErrorContains(t, err, "q2")
	requireBTSCoverageCiphertextPayloadEqual(t, wantInput, input)
	requireBTSCoverageCiphertextPayloadEqual(t, fullOutputBefore, fullOutput)
}

func TestPublicFastBTSLogN16CiphertextMulUsesSharedCoreAndRelinearizeBoundary(t *testing.T) {
	params := btsCoverageLogN16Parameters(t)
	level := params.MaxLevel() - 3
	evaluator := NewEvaluator(params, nil)
	recorder := &recordingMulCore{delegate: fastcore.NewMulWorkspace(params.RingQ(), fastcore.MaxQPrefixWidth)}
	evaluator.mulCore = recorder
	input := btsCoverageCiphertext(t, params, level, 271, false)
	product := NewCiphertext(params, 2, level)
	product.IsNTT = true
	require.NoError(t, evaluator.Mul(input, input, product))
	require.Equal(t, []mulCoreCall{{rows: fastcore.MaxQPrefixWidth, relin: false}}, recorder.calls)
	require.Equal(t, 2, product.Degree())
	require.Equal(t, level, product.Level())
	require.True(t, product.Scale.Equal(input.Scale.Mul(input.Scale)))
	for row := 0; row < fastcore.MaxQPrefixWidth; row++ {
		for _, value := range product.Value[1].Coeffs[row] {
			require.Zero(t, value)
		}
		for _, value := range product.Value[2].Coeffs[row] {
			require.Zero(t, value)
		}
	}

	relinearized := fastcore.NewCompactCiphertext(params, 1, level)
	require.NoError(t, evaluator.Relinearize(product, relinearized))
	require.Equal(t, 1, relinearized.Degree())
	require.Equal(t, level, relinearized.Level())
	require.NoError(t, evaluator.Rescale(relinearized, relinearized))
	require.Equal(t, level-params.LevelsConsumedPerRescaling(), relinearized.Level())
}
