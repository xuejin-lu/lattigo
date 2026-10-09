package fastcore

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func rescaleTransactionTestCiphertext(t *testing.T, ringQ *ring.Ring, level int) *rlwe.Ciphertext {
	t.Helper()
	polys := []ring.Poly{ring.NewPoly(ringQ.N(), level), ring.NewPoly(ringQ.N(), level)}
	ct, err := rlwe.NewCiphertextAtLevelFromPoly(level, polys)
	require.NoError(t, err)
	ct.IsNTT = true
	ct.Scale = rlwe.NewScale(1 << 30)
	ct.MetaData.IsBatched = true
	ct.MetaData.IsBitReversed = true
	ct.MetaData.LogDimensions = ring.Dimensions{Rows: 2, Cols: 3}
	return ct
}

func fillRescaleTransactionInput(t *testing.T, ringQ *ring.Ring, ct *rlwe.Ciphertext) {
	t.Helper()
	divisor := ringQ.SubRings[ct.Level()].Modulus
	for component := range ct.Value {
		value := uint64(component) * divisor
		for row := 0; row <= ct.Level(); row++ {
			modulus := ringQ.SubRings[row].Modulus
			ct.Value[component].Coeffs[row][0] = value % modulus
			ringQ.SubRings[row].NTT(ct.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
		}
	}
}

func requireRescaleCiphertextUnchanged(t *testing.T, before, after *rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, before.Level(), after.Level())
	require.True(t, before.Scale.Equal(after.Scale))
	require.Equal(t, before.MetaData, after.MetaData)
	require.Equal(t, before.IsNTT, after.IsNTT)
	require.Equal(t, before.IsMontgomery, after.IsMontgomery)
	require.Equal(t, before.Value, after.Value)
	for component := range before.Value {
		require.Equal(t, len(before.Value[component].Coeffs), len(after.Value[component].Coeffs))
		for row := range before.Value[component].Coeffs {
			require.Equal(t, cap(before.Value[component].Coeffs[row]), cap(after.Value[component].Coeffs[row]))
		}
	}
}

func TestRescaleCapacityFailureAfterEarlierComponentStagingIsTransactional(t *testing.T) {
	const n, level = 16, 2
	generator := ring.NewNTTFriendlyPrimesGenerator(40, 2*n)
	moduli, err := generator.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	ringQ, err := ring.NewRing(n, moduli)
	require.NoError(t, err)

	for _, inPlace := range []bool{false, true} {
		t.Run(map[bool]string{false: "out_of_place", true: "in_place"}[inPlace], func(t *testing.T) {
			input := rescaleTransactionTestCiphertext(t, ringQ, level)
			fillRescaleTransactionInput(t, ringQ, input)
			inputBefore := input.CopyNew()
			output := rescaleTransactionTestCiphertext(t, ringQ, level)
			for component := range output.Value {
				for row := range output.Value[component].Coeffs {
					for coefficient := range output.Value[component].Coeffs[row] {
						output.Value[component].Coeffs[row][coefficient] = uint64(13 + component + row + coefficient)
					}
				}
			}
			outputBefore := output.CopyNew()
			if inPlace {
				output = input
				outputBefore = inputBefore
			}

			workspace := NewRescaleWorkspace(ringQ, 3)
			// Force the second component's rounded magnitude (1) to fail after
			// the first, zero-valued component has already been staged.
			workspace.half[1] = Uint192{}
			err := workspace.ApplyRows(ringQ, input, output, 1, 3)
			var capacityErr *QPrefixCapacityError
			require.ErrorAs(t, err, &capacityErr)
			require.Equal(t, level-1, capacityErr.Level)
			require.Equal(t, 1, capacityErr.Component)
			requireRescaleCiphertextUnchanged(t, inputBefore, input)
			if !inPlace {
				requireRescaleCiphertextUnchanged(t, outputBefore, output)
			}
			require.Equal(t, uint64(0), workspace.staged[0].Coeffs[0][0])
		})
	}
}

func TestRescaleToNoOpRejectsNoncanonicalInputBeforeOutputMutation(t *testing.T) {
	const n, level = 16, 2
	generator := ring.NewNTTFriendlyPrimesGenerator(40, 2*n)
	moduli, err := generator.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	ringQ, err := ring.NewRing(n, moduli)
	require.NoError(t, err)
	input := rescaleTransactionTestCiphertext(t, ringQ, level)
	fillRescaleTransactionInput(t, ringQ, input)
	input.Scale = rlwe.NewScale(1)
	input.Value[0].Coeffs[0][0] = ringQ.SubRings[0].Modulus
	output := rescaleTransactionTestCiphertext(t, ringQ, level)
	for component := range output.Value {
		for row := range output.Value[component].Coeffs {
			for coefficient := range output.Value[component].Coeffs[row] {
				output.Value[component].Coeffs[row][coefficient] = uint64(17 + component + row + coefficient)
			}
		}
	}
	before := output.CopyNew()

	workspace := NewRescaleWorkspace(ringQ, 3)
	err = workspace.ApplyToRows(ringQ, input, rlwe.NewScale(1<<60), 3, output)
	require.ErrorContains(t, err, "not a canonical residue")
	requireRescaleCiphertextUnchanged(t, before, output)
}

func TestRescaleToNoOpRetainsOnlyExplicitSourceRows(t *testing.T) {
	const n, level = 16, 3
	generator := ring.NewNTTFriendlyPrimesGenerator(40, 2*n)
	moduli, err := generator.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	ringQ, err := ring.NewRing(n, moduli)
	require.NoError(t, err)

	for _, inPlace := range []bool{false, true} {
		t.Run(map[bool]string{false: "out_of_place", true: "in_place"}[inPlace], func(t *testing.T) {
			input := rescaleTransactionTestCiphertext(t, ringQ, level)
			fillRescaleTransactionInput(t, ringQ, input)
			wantRows := make([][]uint64, len(input.Value))
			for component := range input.Value {
				wantRows[component] = append([]uint64(nil), input.Value[component].Coeffs[0]...)
				for row := 1; row <= level; row++ {
					input.Value[component].Coeffs[row] = nil // unauthorized rows must not be read
				}
			}
			output := rescaleTransactionTestCiphertext(t, ringQ, level)
			if inPlace {
				output = input
			}

			workspace := NewRescaleWorkspace(ringQ, 4)
			require.NoError(t, workspace.ApplyToRows(ringQ, input, rlwe.NewScale(1<<40), 1, output))
			require.Equal(t, level, output.Level(), "a no-op must not consume Level")
			require.True(t, rlwe.NewScale(1<<30).Equal(output.Scale))
			for component := range output.Value {
				require.Equal(t, wantRows[component], output.Value[component].Coeffs[0])
				for row := 1; row <= level; row++ {
					require.Nil(t, output.Value[component].Coeffs[row], "unauthorized q%d backing must remain absent", row)
				}
			}
		})
	}
}

func TestRescaleRowsRetainsOnlyAuthorizedTargetPrefix(t *testing.T) {
	const n, level = 16, 3
	generator := ring.NewNTTFriendlyPrimesGenerator(40, 2*n)
	moduli, err := generator.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	ringQ, err := ring.NewRing(n, moduli)
	require.NoError(t, err)
	input := rescaleTransactionTestCiphertext(t, ringQ, level)
	for component := range input.Value {
		for row := 1; row <= level; row++ {
			input.Value[component].Coeffs[row] = nil // only q0 is authorized
		}
	}
	output := rescaleTransactionTestCiphertext(t, ringQ, level-1)
	workspace := NewRescaleWorkspace(ringQ, 4)
	require.NoError(t, workspace.ApplyRows(ringQ, input, output, 1, 1))
	require.Equal(t, level-1, output.Level())
	for component := range output.Value {
		require.NotNil(t, output.Value[component].Coeffs[0])
		require.Nil(t, output.Value[component].Coeffs[1])
		require.Nil(t, output.Value[component].Coeffs[2])
	}
}

func TestRescaleWorkspaceReusesHigherDegreeStaging(t *testing.T) {
	const n, level, degree = 16, 2, 2
	generator := ring.NewNTTFriendlyPrimesGenerator(40, 2*n)
	moduli, err := generator.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	ringQ, err := ring.NewRing(n, moduli)
	require.NoError(t, err)
	inputPolys := make([]ring.Poly, degree+1)
	outputPolys := make([]ring.Poly, degree+1)
	for i := range inputPolys {
		inputPolys[i] = ring.NewPoly(n, level)
		outputPolys[i] = ring.NewPoly(n, level)
	}
	input, err := rlwe.NewCiphertextAtLevelFromPoly(level, inputPolys)
	require.NoError(t, err)
	input.IsNTT = true
	input.Scale = rlwe.NewScale(1 << 30)
	output, err := rlwe.NewCiphertextAtLevelFromPoly(level, outputPolys)
	require.NoError(t, err)
	workspace := NewRescaleWorkspace(ringQ, level+1)
	require.NoError(t, workspace.ApplyRows(ringQ, input, output, 1, level+1))
	require.Len(t, output.Value, degree+1)
	stagedBacking := &workspace.staged[degree].Coeffs[0][0]
	require.NoError(t, workspace.ApplyRows(ringQ, input, output, 1, level+1))
	require.Same(t, stagedBacking, &workspace.staged[degree].Coeffs[0][0], "higher-degree staging backing must be reused")
}
