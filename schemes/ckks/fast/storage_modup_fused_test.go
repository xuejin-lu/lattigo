package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFusedLevel0ModUpToCompactLogicalMatchesPrivateFReference(t *testing.T) {
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
			targetLevel := params.MaxLevel()
			maintained := MaintainedLimbCount(&params, targetLevel)
			for _, inputNTT := range []bool{false, true} {
				for _, targetNTT := range []bool{false, true} {
					t.Run(domainTestName(inputNTT, targetNTT), func(t *testing.T) {
						source := fusedLevel0LogicalSource(t, params, inputNTT)
						before := source.CopyNew()
						want := fusedLevel0PrivateFReference(t, params, source, targetLevel, FastCiphertextDomain{IsNTT: targetNTT})
						got, err := FusedLevel0ModUpToCompactLogical(params, source, targetLevel, FastCiphertextDomain{IsNTT: targetNTT})
						require.NoError(t, err)

						require.Equal(t, targetLevel, got.Level())
						require.Equal(t, source.Degree(), got.Degree())
						require.Equal(t, params.N(), got.N())
						require.True(t, want.Scale.Equal(got.Scale))
						require.Equal(t, targetNTT, got.IsNTT)
						require.False(t, got.IsMontgomery)
						require.Equal(t, want.IsBatched, got.IsBatched)
						require.Equal(t, want.IsBitReversed, got.IsBitReversed)
						require.Equal(t, want.LogDimensions, got.LogDimensions)
						for component := range got.Value {
							for logicalIndex := 0; logicalIndex < maintained; logicalIndex++ {
								require.Equal(t, want.Value[component].Coeffs[logicalIndex], got.Value[component].Coeffs[logicalIndex],
									"component=%d q-index=%d", component, logicalIndex)
							}
							for logicalIndex := maintained; logicalIndex <= targetLevel; logicalIndex++ {
								require.Nil(t, got.Value[component].Coeffs[logicalIndex], "dormant q-index=%d", logicalIndex)
							}
						}

						for component := range source.Value {
							require.Equal(t, before.Value[component].Coeffs[0], source.Value[component].Coeffs[0], "input component=%d was mutated", component)
						}
						require.True(t, before.Scale.Equal(source.Scale))
						require.Equal(t, before.IsNTT, source.IsNTT)
						require.Equal(t, before.IsMontgomery, source.IsMontgomery)
						require.Equal(t, before.IsBatched, source.IsBatched)
						require.Equal(t, before.IsBitReversed, source.IsBitReversed)
						require.Equal(t, before.LogDimensions, source.LogDimensions)
					})
				}
			}
		})
	}
}

func TestFusedLevel0ModUpToCompactLogicalNoFullRNS(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	targetLevel := params.MaxLevel()
	maintained := MaintainedLimbCount(&params, targetLevel)
	source := rlwe.NewCiphertext(&params, 1, 0)
	source.Scale = rlwe.NewScale(1 << 40)
	source.IsNTT = true
	output, err := FusedLevel0ModUpToCompactLogical(params, source, targetLevel, FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	require.Equal(t, targetLevel, output.Level())
	for component := range output.Value {
		for logicalIndex, row := range output.Value[component].Coeffs {
			if logicalIndex < maintained {
				require.Len(t, row, params.N())
			} else {
				require.Nil(t, row, "fused ModUp materialized dormant q-index=%d", logicalIndex)
			}
		}
	}
}

func TestFusedLevel0ModUpToCompactLogicalRejectsInvalidInput(t *testing.T) {
	params := storageBoundaryParameters(t)
	valid := fusedLevel0LogicalSource(t, params, false)
	_, err := FusedLevel0ModUpToCompactLogical(params, nil, params.MaxLevel(), FastCiphertextDomain{})
	require.Error(t, err)
	_, err = FusedLevel0ModUpToCompactLogical(params, valid, 0, FastCiphertextDomain{})
	require.Error(t, err)
	_, err = FusedLevel0ModUpToCompactLogical(params, valid, params.MaxLevel()+1, FastCiphertextDomain{})
	require.Error(t, err)
	_, err = FusedLevel0ModUpToCompactLogical(params, valid, params.MaxLevel(), FastCiphertextDomain{IsMontgomery: true})
	require.Error(t, err)

	montgomery := valid.CopyNew()
	montgomery.IsMontgomery = true
	_, err = FusedLevel0ModUpToCompactLogical(params, montgomery, params.MaxLevel(), FastCiphertextDomain{})
	require.Error(t, err)

	nonCanonical := valid.CopyNew()
	nonCanonical.IsNTT = false
	nonCanonical.Value[0].Coeffs[0][0] = params.Q()[0]
	_, err = FusedLevel0ModUpToCompactLogical(params, nonCanonical, params.MaxLevel(), FastCiphertextDomain{})
	require.Error(t, err)
}

func fusedLevel0LogicalSource(t *testing.T, params ckks.Parameters, inputNTT bool) *rlwe.Ciphertext {
	t.Helper()
	source := rlwe.NewCiphertext(&params, 2, 0)
	source.Scale = rlwe.NewScale(1 << 20)
	source.IsNTT = false
	source.IsBatched = false
	source.IsBitReversed = true
	source.LogDimensions = ring.Dimensions{Rows: 2, Cols: 3}
	q0 := params.Q()[0]
	half := q0 >> 1
	values := []uint64{0, 1, half - 1, half, half + 1, q0 - 1, q0 - 2, half - 7, half + 8}
	for component := range source.Value {
		row := source.Value[component].Coeffs[0]
		for coefficient := range row {
			row[coefficient] = values[(coefficient+component)%len(values)]
		}
		if inputNTT {
			params.RingQ().SubRings[0].NTT(row, row)
		}
	}
	source.IsNTT = inputNTT
	return source
}

func fusedLevel0PrivateFReference(t *testing.T, params ckks.Parameters, source *rlwe.Ciphertext, targetLevel int, target FastCiphertextDomain) *rlwe.Ciphertext {
	t.Helper()
	private, err := ImportLevel0(params, source, 3, FastCiphertextDomain{IsNTT: source.IsNTT})
	require.NoError(t, err)
	private, err = FastStorageModUpLevel0(private, targetLevel)
	require.NoError(t, err)
	output, err := private.ExportToCompactLogical(target)
	require.NoError(t, err)
	return output
}

func domainTestName(inputNTT, targetNTT bool) string {
	inputDomain, outputDomain := "coeff", "coeff"
	if inputNTT {
		inputDomain = "ntt"
	}
	if targetNTT {
		outputDomain = "ntt"
	}
	return inputDomain + "-to-" + outputDomain
}
