package ckks

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

func plaintextAddTestParameters(t *testing.T) Parameters {
	t.Helper()
	logQ := []int{42}
	for i := 0; i < 15; i++ {
		logQ = append(logQ, 32)
	}
	params, err := NewParametersFromLiteral(ParametersLiteral{LogN: 10, LogQ: logQ, LogDefaultScale: 32})
	require.NoError(t, err)
	return params
}

func plaintextAddTestInputs(t *testing.T, params Parameters, ctLevel, ptLevel int) (*rlwe.Ciphertext, *rlwe.Plaintext, []complex128) {
	t.Helper()
	values := fastAddSubTestValues(1)
	bias := fastAddSubTestValues(5)
	encoder := NewEncoder(params)
	message := NewPlaintext(params, ctLevel)
	message.LogDimensions = ring.Dimensions{Cols: fastRotateTestLogSlots}
	// Match stem.square's non-power-of-two scale after L14 -> L13 Rescale.
	message.Scale = rlwe.NewScale(math.Exp2(64)).Div(rlwe.NewScale(params.Q()[14]))
	require.NoError(t, encoder.Encode(values, message))
	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	full, err := rlwe.NewEncryptor(params, sk).EncryptNew(message)
	require.NoError(t, err)
	ct := fastcore.NewCompactCiphertext(params, 1, ctLevel)
	*ct.MetaData = *full.MetaData
	rows, err := fastcore.QPrefixWidth(ctLevel)
	require.NoError(t, err)
	for component := range ct.Value {
		fastcore.CopyPrefixRowsUnchecked(rows, full.Value[component], ct.Value[component])
	}
	pt := NewPlaintext(params, ptLevel)
	pt.LogDimensions = message.LogDimensions
	pt.Scale = ct.Scale
	require.NoError(t, encoder.Encode(bias, pt))
	want := make([]complex128, len(values))
	for i := range want {
		want[i] = values[i] + bias[i]
	}
	return ct, pt, want
}

func TestPublicFastPlaintextAddUsesQPrefixAndDecodes(t *testing.T) {
	params := plaintextAddTestParameters(t)
	for _, levels := range [][3]int{{13, 13, 13}, {13, 15, 13}, {15, 13, 15}, {13, 13, 8}, {3, 2, 3}, {2, 3, 2}, {1, 1, 1}} {
		for _, inPlace := range []bool{false, true} {
			if inPlace && levels[2] != levels[0] {
				continue
			}
			t.Run(fmt.Sprintf("ct%d-pt%d-out%d/inplace=%t", levels[0], levels[1], levels[2], inPlace), func(t *testing.T) {
				ct, pt, want := plaintextAddTestInputs(t, params, levels[0], levels[1])
				ctBefore := ct.CopyNew()
				ptBefore := pt.CopyNew()
				out := NewCiphertext(params, 1, levels[2])
				if inPlace {
					out = ct
				}
				eval := NewEvaluator(params, nil)
				recorder := &recordingAddSubCore{delegate: eval.addSubCore}
				eval.addSubCore = recorder
				require.NoError(t, eval.Add(ct, pt, out))
				level := min(levels[0], levels[1], levels[2])
				rows, err := fastcore.QPrefixWidth(level)
				require.NoError(t, err)
				require.Equal(t, []addSubCoreCall{{rows: rows}}, recorder.calls)
				require.Equal(t, level, out.Level())
				require.Equal(t, 1, out.Degree())
				require.True(t, ctBefore.Scale.Equal(out.Scale))
				require.Equal(t, ctBefore.IsNTT, out.IsNTT)
				require.Equal(t, ctBefore.IsMontgomery, out.IsMontgomery)
				for row := 0; row < rows; row++ {
					q := params.Q()[row]
					for i, value := range ctBefore.Value[0].Coeffs[row] {
						require.Equal(t, (value+ptBefore.Value.Coeffs[row][i])%q, out.Value[0].Coeffs[row][i])
					}
					require.Equal(t, ctBefore.Value[1].Coeffs[row], out.Value[1].Coeffs[row])
				}
				for component := range out.Value {
					for row := rows; row <= level; row++ {
						require.Nil(t, out.Value[component].Coeffs[row])
					}
				}
				require.True(t, pt.Equal(ptBefore), "plaintext must not be mutated")
				if !inPlace {
					require.True(t, ct.Equal(ctBefore), "out-of-place Add must not mutate its ciphertext input")
				}
				// Fast output has zero c1; ordinary DecryptNew/Decode provides the oracle boundary.
				sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
				decoded := decodeFastAddSubTestValues(t, params, NewEncoder(params), rlwe.NewDecryptor(params, sk), out)
				var maxError float64
				for i := range want {
					maxError = math.Max(maxError, cmplx.Abs(decoded[i]-want[i]))
				}
				t.Logf("max complex error %.12g (limit 1e-6)", maxError)
				require.LessOrEqual(t, maxError, 1e-6)
			})
		}
	}
}

func TestPublicFastPlaintextAddIgnoresDormantRows(t *testing.T) {
	params := plaintextAddTestParameters(t)
	for _, missing := range []bool{false, true} {
		t.Run(fmt.Sprintf("missing=%t", missing), func(t *testing.T) {
			ct, pt, want := plaintextAddTestInputs(t, params, 13, 13)
			for row := fastcore.MaxQPrefixWidth; row <= pt.Level(); row++ {
				if missing {
					pt.Value.Coeffs[row] = nil
				} else {
					for i := range pt.Value.Coeffs[row] {
						pt.Value.Coeffs[row][i] = ^uint64(0)
					}
				}
			}
			before := pt.CopyNew()
			eval := NewEvaluator(params, nil)
			recorder := &recordingAddSubCore{delegate: eval.addSubCore}
			eval.addSubCore = recorder
			require.NoError(t, eval.Add(ct, pt, ct))
			require.Equal(t, []addSubCoreCall{{rows: 4}}, recorder.calls)
			require.True(t, pt.Equal(before))
			sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
			requireFastChainOracle(t, params, NewEncoder(params), rlwe.NewDecryptor(params, sk), ct, want, 1e-6, "poisoned/dormant plaintext rows")
		})
	}
}

func TestPublicFastPlaintextAddRejectsInvalidInputsWithoutMutation(t *testing.T) {
	params := plaintextAddTestParameters(t)
	for _, name := range []string{"nil-pt", "nil-ct", "nil-output", "pt-metadata", "ct-metadata", "output-metadata", "empty-pt", "empty-ct", "empty-output", "scale", "NTT", "Montgomery", "batching", "pt-row", "ct-row", "output-row", "nonzero-c1", "residue", "ring-degree"} {
		t.Run(name, func(t *testing.T) {
			ct, pt, _ := plaintextAddTestInputs(t, params, 13, 13)
			out := fastcore.NewCompactCiphertext(params, 1, 13)
			switch name {
			case "nil-pt":
				pt = nil
			case "nil-ct":
				ct = nil
			case "nil-output":
				out = nil
			case "pt-metadata":
				pt.MetaData = nil
			case "ct-metadata":
				ct.MetaData = nil
			case "output-metadata":
				out.MetaData = nil
			case "empty-pt":
				pt.Element.Value = nil
			case "empty-ct":
				ct.Value = nil
			case "empty-output":
				out.Value = nil
			case "scale":
				pt.Scale = rlwe.NewScale(1)
			case "NTT":
				pt.IsNTT = false
			case "Montgomery":
				pt.IsMontgomery = true
			case "batching":
				pt.IsBatched = !ct.IsBatched
			case "pt-row":
				pt.Value.Coeffs[2] = nil
			case "ct-row":
				ct.Value[0].Coeffs[2] = nil
			case "output-row":
				out.Value[0].Coeffs[2] = nil
			case "nonzero-c1":
				ct.Value[1].Coeffs[0][0] = 1
			case "residue":
				pt.Value.Coeffs[0][0] = params.Q()[0]
			case "ring-degree":
				pt.Value.Coeffs[0] = pt.Value.Coeffs[0][:params.N()/2]
			}
			// Snapshot only existing prefix rows, retaining nil slice identity.
			var snapshots [][][]uint64
			var polys []*ring.Poly
			for _, obj := range []*rlwe.Ciphertext{ct, out} {
				if obj != nil {
					for i := range obj.Value {
						polys = append(polys, &obj.Value[i])
					}
				}
			}
			if pt != nil && len(pt.Element.Value) == 1 {
				polys = append(polys, &pt.Element.Value[0])
			}
			for _, poly := range polys {
				rows := make([][]uint64, min(4, len(poly.Coeffs)))
				for row := range rows {
					if poly.Coeffs[row] != nil {
						rows[row] = append([]uint64{}, poly.Coeffs[row]...)
					}
				}
				snapshots = append(snapshots, rows)
			}
			eval := NewEvaluator(params, nil)
			require.NotPanics(t, func() { require.Error(t, eval.Add(ct, pt, out)) })
			for i, poly := range polys {
				for row := range snapshots[i] {
					require.Equal(t, snapshots[i][row], poly.Coeffs[row])
				}
			}
		})
	}
}

// Full native ciphertexts outside the zero-secret contract retain the old
// generic Add behavior; this fork-local test is not a genuine Standard baseline.
func TestPublicFastPlaintextAddPreservesNativeFullCiphertextFallback(t *testing.T) {
	params := plaintextAddTestParameters(t)
	native := nativeEncryptionParameters{params}
	kgen := rlwe.NewKeyGenerator(native)
	sk := kgen.GenSecretKeyNew()
	encoder := NewEncoder(params)
	message := NewPlaintext(params, 13)
	message.LogDimensions = ring.Dimensions{Cols: fastRotateTestLogSlots}
	values, bias := fastAddSubTestValues(1), fastAddSubTestValues(5)
	require.NoError(t, encoder.Encode(values, message))
	ct, err := rlwe.NewEncryptor(native, sk).EncryptNew(message)
	require.NoError(t, err)
	require.False(t, fastQPolyIsZero(ct.Value[1], 4))
	pt := NewPlaintext(params, 13)
	pt.LogDimensions, pt.Scale = message.LogDimensions, message.Scale
	require.NoError(t, encoder.Encode(bias, pt))
	eval := NewEvaluator(params, nil)
	recorder := &recordingAddSubCore{delegate: eval.addSubCore}
	eval.addSubCore = recorder
	require.NoError(t, eval.Add(ct, pt, ct))
	require.Empty(t, recorder.calls)
	decoded := decodeFastAddSubTestValues(t, params, encoder, rlwe.NewDecryptor(native, sk), ct)
	want := make([]complex128, len(values))
	for i := range want {
		want[i] = values[i] + bias[i]
	}
	requireComplexValuesNear(t, decoded, want, 1e-6)
}
