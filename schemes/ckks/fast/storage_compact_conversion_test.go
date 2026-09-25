package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFastStorageExportToCompactLogicalStructureMetadataAndDomains(t *testing.T) {
	profiles := []struct {
		name   string
		params ckks.Parameters
	}{
		{name: "q01", params: storageBoundaryParameters(t)},
		{name: "q012", params: compactLogicalQ012Parameters(t)},
	}
	for _, profile := range profiles {
		t.Run(profile.name, func(t *testing.T) {
			params := profile.params
			level := params.MaxLevel()
			maintained := MaintainedLimbCount(&params, level)
			require.Greater(t, level+1, maintained, "fixture must have dormant logical rows")
			values := compactLogicalQFixtureValues(params)
			for _, inputNTT := range []bool{false, true} {
				input := newFastStorageCiphertextFromBig(t, params, values, level, 3, inputNTT, nil)
				input.metadata.Scale = rlwe.NewScale(123456789)
				input.metadata.IsBatched = false
				input.metadata.IsBitReversed = true
				input.metadata.LogDimensions = ring.Dimensions{Rows: 2, Cols: 3}
				inputMetadata := input.MetaData()
				for _, targetNTT := range []bool{false, true} {
					output, err := input.ExportToCompactLogical(FastCiphertextDomain{IsNTT: targetNTT})
					require.NoError(t, err)
					require.Equal(t, level, output.Level())
					require.Equal(t, input.Degree(), output.Degree())
					require.True(t, input.Scale().Equal(output.Scale))
					require.Equal(t, targetNTT, output.IsNTT)
					require.False(t, output.IsMontgomery)
					require.Equal(t, inputMetadata.IsBatched, output.IsBatched)
					require.Equal(t, inputMetadata.IsBitReversed, output.IsBitReversed)
					require.Equal(t, inputMetadata.LogDimensions, output.LogDimensions)
					for component := range output.Value {
						require.Len(t, output.Value[component].Coeffs, level+1)
						for logicalIndex, row := range output.Value[component].Coeffs {
							if logicalIndex < maintained {
								require.Len(t, row, params.N(), "maintained logical row %d", logicalIndex)
							} else {
								require.Nil(t, row, "dormant logical row %d must remain unmaterialized", logicalIndex)
							}
						}
					}
				}
			}
		})
	}
}

func TestFastStorageExportToCompactLogicalExactness(t *testing.T) {
	profiles := []struct {
		name   string
		params ckks.Parameters
	}{
		{name: "q01", params: storageBoundaryParameters(t)},
		{name: "q012", params: compactLogicalQ012Parameters(t)},
	}
	for _, profile := range profiles {
		t.Run(profile.name, func(t *testing.T) {
			params := profile.params
			level := params.MaxLevel()
			values := compactLogicalQFixtureValues(params)
			maintained := MaintainedLimbCount(&params, level)
			want := compactLogicalQExpected(values, params, maintained)
			for _, inputNTT := range []bool{false, true} {
				for _, targetNTT := range []bool{false, true} {
					input := newFastStorageCiphertextFromBig(t, params, values, level, 3, inputNTT, nil)
					output, err := input.ExportToCompactLogical(FastCiphertextDomain{IsNTT: targetNTT})
					require.NoError(t, err)
					for component := range output.Value {
						for logicalIndex := 0; logicalIndex < maintained; logicalIndex++ {
							row := output.Value[component].Coeffs[logicalIndex]
							if targetNTT {
								coefficients := make([]uint64, params.N())
								params.RingQ().SubRings[logicalIndex].INTT(row, coefficients)
								row = coefficients
							}
							require.Equal(t, want[component][logicalIndex], row,
								"inputNTT=%t targetNTT=%t component=%d q-index=%d", inputNTT, targetNTT, component, logicalIndex)
						}
					}
				}
			}
		})
	}
}

func TestFastStorageExportToCompactLogicalDoesNotMaterializeFullRNS(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	level := params.MaxLevel()
	maintained := MaintainedLimbCount(&params, level)
	require.Greater(t, level+1, maintained)
	values := compactLogicalQFixtureValues(params)
	input := newFastStorageCiphertextFromBig(t, params, values, level, 3, true, nil)

	output, err := input.ExportToCompactLogical(FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	require.Equal(t, level, output.Level())
	for component := range output.Value {
		for logicalIndex, row := range output.Value[component].Coeffs {
			if logicalIndex < maintained {
				require.Len(t, row, params.N())
			} else {
				require.Nil(t, row, "full-RNS logical row %d was unexpectedly materialized", logicalIndex)
			}
		}
	}
}

func compactLogicalQ012Parameters(t *testing.T) ckks.Parameters {
	t.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(45, 2*2*16)
	q3, err := generator.NextAlternatingPrimes(1)
	require.NoError(t, err)
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		Q:               []uint64{72057594037616641, 549755731969, 549756026881, q3[0]},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	require.True(t, q012Enabled(&params), "fixture must exercise the maintained q012 policy")
	return params
}

func compactLogicalQFixtureValues(params ckks.Parameters) [][]*big.Int {
	q0 := new(big.Int).SetUint64(params.Q()[0])
	return [][]*big.Int{
		{
			new(big.Int).Add(big.NewInt(7), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(13))),
			new(big.Int).Sub(big.NewInt(-11), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(5))),
			new(big.Int).Add(new(big.Int).Rsh(new(big.Int).Set(q0), 1), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(3))),
		},
		{
			new(big.Int).Neg(new(big.Int).Add(big.NewInt(23), new(big.Int).Mul(new(big.Int).Set(q0), big.NewInt(17)))),
			big.NewInt(0),
		},
	}
}

func compactLogicalQExpected(values [][]*big.Int, params ckks.Parameters, maintained int) [][][]uint64 {
	want := make([][][]uint64, len(values))
	for component := range values {
		want[component] = make([][]uint64, maintained)
		for logicalIndex := 0; logicalIndex < maintained; logicalIndex++ {
			q := new(big.Int).SetUint64(params.Q()[logicalIndex])
			want[component][logicalIndex] = make([]uint64, params.N())
			for coefficient, value := range values[component] {
				want[component][logicalIndex][coefficient] = new(big.Int).Mod(new(big.Int).Set(value), q).Uint64()
			}
		}
	}
	return want
}
