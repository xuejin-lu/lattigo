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

func TestFastRescaleStagesHigherDegreeComponents(t *testing.T) {
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

	require.NoError(t, fastEval.RescaleQPrefixRows(fastInput, 4, fastOutput))
	require.Len(t, fastOutput.Value, 3)
}
