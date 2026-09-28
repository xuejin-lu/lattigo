package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type fastRescalePolyShape struct {
	rowCount      int
	rowCapacity   int
	rowLengths    []int
	rowCapacities []int
}

type fastRescaleCiphertextShape struct {
	valueCount    int
	valueCapacity int
	polys         []fastRescalePolyShape
}

type fastRescaleSnapshot struct {
	value        *rlwe.Ciphertext
	shape        fastRescaleCiphertextShape
	level        int
	scale        rlwe.Scale
	metadata     rlwe.MetaData
	isNTT        bool
	isMontgomery bool
}

func snapshotFastRescaleCiphertext(ct *rlwe.Ciphertext) fastRescaleSnapshot {
	shape := fastRescaleCiphertextShape{valueCount: len(ct.Value), valueCapacity: cap(ct.Value), polys: make([]fastRescalePolyShape, len(ct.Value))}
	for i := range ct.Value {
		poly := ct.Value[i]
		rowShape := fastRescalePolyShape{rowCount: len(poly.Coeffs), rowCapacity: cap(poly.Coeffs), rowLengths: make([]int, len(poly.Coeffs)), rowCapacities: make([]int, len(poly.Coeffs))}
		for row := range poly.Coeffs {
			rowShape.rowLengths[row] = len(poly.Coeffs[row])
			rowShape.rowCapacities[row] = cap(poly.Coeffs[row])
		}
		shape.polys[i] = rowShape
	}
	return fastRescaleSnapshot{
		value: ct.CopyNew(), shape: shape, level: ct.Level(), scale: ct.Scale,
		metadata: *ct.MetaData, isNTT: ct.IsNTT, isMontgomery: ct.IsMontgomery,
	}
}

func requireFastRescaleUnchanged(t *testing.T, before fastRescaleSnapshot, after *rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, before.level, after.Level(), "logical Level changed")
	require.True(t, before.scale.Equal(after.Scale), "Scale changed")
	require.Equal(t, before.metadata, *after.MetaData, "metadata changed")
	require.Equal(t, before.isNTT, after.IsNTT, "NTT flag changed")
	require.Equal(t, before.isMontgomery, after.IsMontgomery, "Montgomery flag changed")
	require.Equal(t, before.value.Value, after.Value, "coefficient backing changed")
	require.Equal(t, before.shape, fastRescaleShape(after), "ciphertext/storage shape or backing capacities changed")
}

func fastRescaleShape(ct *rlwe.Ciphertext) fastRescaleCiphertextShape {
	shape := fastRescaleCiphertextShape{valueCount: len(ct.Value), valueCapacity: cap(ct.Value), polys: make([]fastRescalePolyShape, len(ct.Value))}
	for i := range ct.Value {
		poly := ct.Value[i]
		rowShape := fastRescalePolyShape{rowCount: len(poly.Coeffs), rowCapacity: cap(poly.Coeffs), rowLengths: make([]int, len(poly.Coeffs)), rowCapacities: make([]int, len(poly.Coeffs))}
		for row := range poly.Coeffs {
			rowShape.rowLengths[row] = len(poly.Coeffs[row])
			rowShape.rowCapacities[row] = cap(poly.Coeffs[row])
		}
		shape.polys[i] = rowShape
	}
	return shape
}

func makeTransactionalRescaleInput(params ckks.Parameters) *rlwe.Ciphertext {
	ct := NewCiphertext(params, 1, 3)
	ct.IsNTT = true
	ct.IsMontgomery = true
	ct.Scale = rlwe.NewScale(1 << 40)
	ct.MetaData.IsBatched = true
	ct.MetaData.IsBitReversed = true
	ct.MetaData.LogDimensions = ring.Dimensions{Rows: 2, Cols: 3}
	q3 := params.Q()[3]
	for component := range ct.Value {
		value := uint64(component+1) * q3
		for row := 0; row < 4; row++ {
			q := params.Q()[row]
			ct.Value[component].Coeffs[row][0] = value % q
			ntt := make([]uint64, params.N())
			params.RingQ().SubRings[row].NTT(ct.Value[component].Coeffs[row], ntt)
			params.RingQ().SubRings[row].MForm(ntt, ntt)
			copy(ct.Value[component].Coeffs[row], ntt)
		}
	}
	return ct
}

func TestFastRescaleCapacityFailureAfterEarlierComponentStagingIsTransactional(t *testing.T) {
	params := q012TestParameters(t)
	for _, inPlace := range []bool{true, false} {
		name := "out-of-place"
		if inPlace {
			name = "in-place"
		}
		t.Run(name, func(t *testing.T) {
			eval := NewEvaluator(params)
			input := makeTransactionalRescaleInput(params)
			output := NewCiphertext(params, 2, 1)
			output.IsNTT = false
			output.IsMontgomery = false
			output.Scale = rlwe.NewScale(12345)
			output.MetaData.IsBatched = true
			output.MetaData.IsBitReversed = true
			output.MetaData.LogDimensions = ring.Dimensions{Rows: 1, Cols: 2}
			for component := range output.Value {
				for row := range output.Value[component].Coeffs {
					for coefficient := range output.Value[component].Coeffs[row] {
						output.Value[component].Coeffs[row][coefficient] = uint64(17 + 13*component + 5*row + coefficient)
					}
				}
			}

			if inPlace {
				output = input
			}
			inputBefore := snapshotFastRescaleCiphertext(input)
			var outputBefore fastRescaleSnapshot
			if !inPlace {
				outputBefore = snapshotFastRescaleCiphertext(output)
			}

			// A valid Rescale over a centered source prefix normally contracts
			// within capacity. Tighten only this evaluator's target-half scratch
			// to inject the later-component capacity failure deterministically:
			// c0 rounds to 1 and stages successfully; c1 rounds to 2 and fails.
			eval.rescaleScratch.half[2] = uint192{lo: 1}
			eval.rescaleScratch.staged[0].Coeffs[0][0] = params.Q()[0] - 1
			err := eval.RescaleQPrefixRows(input, 4, output)
			var capacityErr *QPrefixCapacityError
			require.ErrorAs(t, err, &capacityErr)
			require.Equal(t, 1, capacityErr.Component, "capacity error: %+v", capacityErr)
			require.Equal(t, uint64(1), eval.rescaleScratch.staged[0].Coeffs[0][0], "first component must have staged its computed result")
			requireFastRescaleUnchanged(t, inputBefore, input)
			if !inPlace {
				requireFastRescaleUnchanged(t, outputBefore, output)
			}
		})
	}
}

func TestFastRescaleStagesHigherDegreeComponentsAndReusesThem(t *testing.T) {
	params := q012TestParameters(t)
	level := 3
	fastInput := NewCiphertext(params, 2, level)
	standardInput := ckks.NewCiphertext(params, 2, level)
	fastInput.IsNTT, standardInput.IsNTT = true, true
	fastInput.Scale, standardInput.Scale = rlwe.NewScale(1<<60), rlwe.NewScale(1<<60)
	q3 := params.Q()[3]
	for component := range fastInput.Value {
		for row := 0; row <= level; row++ {
			q := params.Q()[row]
			for coefficient := 0; coefficient < params.N(); coefficient++ {
				value := int64(q3)*int64(component+1) + int64(coefficient%5-2)
				if coefficient%2 == 0 {
					value = -value
				}
				residue := uint64(0)
				if value < 0 {
					magnitude := uint64(-value) % q
					if magnitude != 0 {
						residue = q - magnitude
					}
				} else {
					residue = uint64(value) % q
				}
				fastInput.Value[component].Coeffs[row][coefficient] = residue
				standardInput.Value[component].Coeffs[row][coefficient] = residue
			}
			ntt := make([]uint64, params.N())
			params.RingQ().SubRings[row].NTT(fastInput.Value[component].Coeffs[row], ntt)
			copy(fastInput.Value[component].Coeffs[row], ntt)
			params.RingQ().SubRings[row].NTT(standardInput.Value[component].Coeffs[row], ntt)
			copy(standardInput.Value[component].Coeffs[row], ntt)
		}
	}

	fastEval := NewEvaluator(params)
	fastOutput := NewCiphertext(params, fastInput.Degree(), level-1)
	standardOutput := ckks.NewCiphertext(params, standardInput.Degree(), level-1)
	require.NoError(t, fastEval.RescaleQPrefixRows(fastInput, 4, fastOutput))
	standardEval := ckks.NewEvaluator(params, nil)
	require.NoError(t, standardEval.Rescale(standardInput, standardOutput))
	require.Equal(t, standardOutput.Level(), fastOutput.Level())
	require.True(t, standardOutput.Scale.Equal(fastOutput.Scale))
	require.Len(t, fastOutput.Value, 3)
	for component := range fastOutput.Value {
		for row := 0; row <= fastOutput.Level(); row++ {
			require.Equal(t, standardOutput.Value[component].Coeffs[row], fastOutput.Value[component].Coeffs[row], "component=%d q%d", component, row)
		}
	}

	stagedRow := fastEval.rescaleScratch.staged[2].Coeffs[0]
	require.NoError(t, fastEval.RescaleQPrefixRows(fastInput, 4, fastOutput))
	require.True(t, &stagedRow[0] == &fastEval.rescaleScratch.staged[2].Coeffs[0][0], "higher-degree staging backing must be reused")
}
