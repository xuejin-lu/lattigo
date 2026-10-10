//go:build fastdiag

package fastcore

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/internal/fastdiag"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func TestRescaleApplyRowsFastDiagEmitsActualMaterializationTree(t *testing.T) {
	for sourceRows := 1; sourceRows <= 4; sourceRows++ {
		level := sourceRows
		generator := ring.NewNTTFriendlyPrimesGenerator(40, 2*16)
		moduli, err := generator.NextAlternatingPrimes(level + 1)
		require.NoError(t, err)
		ringQ, err := ring.NewRing(16, moduli)
		require.NoError(t, err)

		for _, montgomery := range []bool{false, true} {
			for _, inPlace := range []bool{false, true} {
				name := fmt.Sprintf("rows_%d/montgomery_%t/in_place_%t", sourceRows, montgomery, inPlace)
				t.Run(name, func(t *testing.T) {
					input := rescaleTransactionTestCiphertext(t, ringQ, level)
					fillRescaleTransactionInput(t, ringQ, input)
					if montgomery {
						for component := range input.Value {
							for row := 0; row < sourceRows; row++ {
								ringQ.SubRings[row].MForm(input.Value[component].Coeffs[row], input.Value[component].Coeffs[row])
							}
						}
						input.IsMontgomery = true
					}

					run := func(trace string) (*rlwe.Ciphertext, []fastdiag.Event) {
						require.NoError(t, fastdiag.ConfigureCSV(trace))
						fastdiag.Reset()
						in := input.CopyNew()
						out := rescaleTransactionTestCiphertext(t, ringQ, level)
						if inPlace {
							out = in
						}
						workspace := NewRescaleWorkspace(ringQ, sourceRows)
						require.NoError(t, workspace.ApplyRows(ringQ, in, out, 1, sourceRows))
						return out.CopyNew(), fastdiag.Events()
					}

					withoutTrace, noTraceEvents := run("")
					require.Empty(t, noTraceEvents)
					withTrace, events := run("rescale")
					require.Equal(t, withoutTrace, withTrace, "diagnostic spans must not alter Rescale output")
					require.Equal(t, sourceRows-1, withTrace.Level())
					require.True(t, withTrace.IsNTT)
					require.Equal(t, montgomery, withTrace.IsMontgomery)
					fastDiagRequireRescaleMaterializationTree(t, events, sourceRows, inPlace)
				})
			}
		}
	}
	require.NoError(t, fastdiag.ConfigureCSV(""))
}

func TestRescaleApplyRowsFastDiagDoesNotFabricateMaterializationOnCapacityFailure(t *testing.T) {
	const n, level = 16, 2
	generator := ring.NewNTTFriendlyPrimesGenerator(40, 2*n)
	moduli, err := generator.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	ringQ, err := ring.NewRing(n, moduli)
	require.NoError(t, err)
	input := rescaleTransactionTestCiphertext(t, ringQ, level)
	fillRescaleTransactionInput(t, ringQ, input)
	output := rescaleTransactionTestCiphertext(t, ringQ, level)
	outputBefore := output.CopyNew()
	workspace := NewRescaleWorkspace(ringQ, 3)
	workspace.half[1] = Uint192{}
	require.NoError(t, fastdiag.ConfigureCSV("rescale"))
	fastdiag.Reset()
	err = workspace.ApplyRows(ringQ, input, output, 1, 3)
	var capacityErr *QPrefixCapacityError
	require.ErrorAs(t, err, &capacityErr)
	requireRescaleCiphertextUnchanged(t, outputBefore, output)

	for _, event := range fastdiag.Events() {
		require.NotEqual(t, "materialization", event.Name, "failed preflight must not emit a completed materialization span")
		require.NotEqual(t, "ntt_montgomery_restore", event.Name, "failed preflight must not emit a restore span")
	}
	require.NoError(t, fastdiag.ConfigureCSV(""))
}

func fastDiagRequireRescaleMaterializationTree(t *testing.T, events []fastdiag.Event, sourceRows int, inPlace bool) {
	t.Helper()
	sequences := make(map[uint64]fastdiag.Event, len(events))
	children := make(map[uint64][]fastdiag.Event, len(events))
	for _, event := range events {
		require.NotZero(t, event.Sequence)
		require.GreaterOrEqual(t, event.ElapsedNS, int64(0))
		require.NotContains(t, sequences, event.Sequence, "duplicate event sequence")
		sequences[event.Sequence] = event
		children[event.ParentSequence] = append(children[event.ParentSequence], event)
	}

	var rescale fastdiag.Event
	parents := 0
	for _, event := range events {
		if event.Scope == fastdiag.Rescale && event.Name == "rescale" {
			parents++
			rescale = event
		}
	}
	require.Equal(t, 1, parents)
	require.NotZero(t, rescale.Sequence)
	for _, event := range events {
		if event.ParentSequence != 0 {
			_, exists := sequences[event.ParentSequence]
			require.True(t, exists, "event %d has orphan parent %d", event.Sequence, event.ParentSequence)
		}
	}

	passes := make(map[string]fastdiag.Event, 2)
	for _, child := range children[rescale.Sequence] {
		require.Equal(t, fastdiag.Rescale, child.Scope)
		passes[child.Name] = child
	}
	require.Len(t, passes, 2)
	preflight, hasPreflight := passes["preflight"]
	require.True(t, hasPreflight)
	materialization, hasMaterialization := passes["materialization"]
	require.True(t, hasMaterialization)
	require.NotNil(t, preflight.Count)
	require.Equal(t, 2, *preflight.Count)
	require.NotNil(t, materialization.Count)
	require.Equal(t, 2, *materialization.Count)
	fastDiagRequireMeasuredClosure(t, rescale, children[rescale.Sequence])
	fastDiagRequireMeasuredClosure(t, preflight, children[preflight.Sequence])
	fastDiagRequireMeasuredClosure(t, materialization, children[materialization.Sequence])

	restores := children[materialization.Sequence]
	require.Len(t, restores, 2)
	restoreComponents := make(map[string]bool, len(restores))
	for _, restore := range restores {
		require.Equal(t, fastdiag.Rescale, restore.Scope)
		require.Equal(t, "ntt_montgomery_restore", restore.Name)
		require.NotEmpty(t, restore.Component)
		require.False(t, restoreComponents[restore.Component], "duplicate restore event for %s", restore.Component)
		restoreComponents[restore.Component] = true
		require.NotNil(t, restore.RowsIn)
		require.Equal(t, sourceRows, *restore.RowsIn)
		require.NotNil(t, restore.RowsOut)
		require.NotNil(t, restore.InPlace)
		require.Equal(t, inPlace, *restore.InPlace)
	}
	require.True(t, restoreComponents["c0"])
	require.True(t, restoreComponents["c1"])

	preflightComponents := make(map[string]map[string]int)
	for _, child := range children[preflight.Sequence] {
		require.Equal(t, fastdiag.Rescale, child.Scope)
		if preflightComponents[child.Component] == nil {
			preflightComponents[child.Component] = make(map[string]int)
		}
		preflightComponents[child.Component][child.Name]++
	}
	require.Len(t, preflightComponents, 2)
	for _, component := range []string{"c0", "c1"} {
		require.Equal(t, map[string]int{"prefix_to_coefficient": 1, "coefficient_loop": 1}, preflightComponents[component])
	}
}

func fastDiagRequireMeasuredClosure(t *testing.T, parent fastdiag.Event, children []fastdiag.Event) {
	t.Helper()
	var childNS int64
	for _, child := range children {
		childNS += child.ElapsedNS
	}
	require.GreaterOrEqual(t, parent.ElapsedNS, childNS,
		"measured %s/%s parent must cover its sequential direct-child intervals; preserve the unnormalized residual", parent.Scope, parent.Name)
}
