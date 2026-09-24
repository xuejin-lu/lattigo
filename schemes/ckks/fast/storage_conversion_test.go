package fast

import (
	"math/big"
	"math/bits"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func storageBoundaryParameters(t *testing.T) ckks.Parameters {
	t.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(45, 2*16)
	q, err := generator.NextAlternatingPrimes(6)
	require.NoError(t, err)
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		Q:               q,
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	require.Equal(t, ring.Standard, params.RingType())
	require.GreaterOrEqual(t, params.MaxLevel(), 5)
	return params
}

func TestFastCiphertextLogicalLevelAndStorageWidthAreIndependent(t *testing.T) {
	params := storageBoundaryParameters(t)
	var uninitialized FastCiphertext
	require.Error(t, uninitialized.ResizeDegree(1))
	for _, test := range []struct {
		level int
		width int
	}{{5, 3}, {1, 3}, {5, 2}} {
		ct, err := NewFastCiphertext(params, 1, test.level, test.width)
		require.NoError(t, err)
		require.Equal(t, test.level, ct.LogicalLevel())
		require.Equal(t, test.width, ct.StorageWidth())
		require.Equal(t, params.N(), ct.N())
		require.Equal(t, 1, ct.Degree())
		for component := range ct.value {
			require.Len(t, ct.value[component].Coeffs, test.width)
			for _, row := range ct.value[component].Coeffs {
				require.Len(t, row, params.N())
			}
		}
	}

	ct, err := NewFastCiphertext(params, 1, 5, 3)
	require.NoError(t, err)
	for component := range ct.value {
		for row := range ct.value[component].Coeffs {
			for k := range ct.value[component].Coeffs[row] {
				ct.value[component].Coeffs[row][k] = uint64(component+row+k+1) % fastStoragePrimes()[row]
			}
		}
	}
	before := cloneStorageRows(ct.value)
	require.NoError(t, ct.SetLogicalLevel(1))
	require.Equal(t, before, ct.value, "logical level changes must not mutate or resize F rows")
	require.Equal(t, 3, ct.StorageWidth())

	copyCT := ct.CopyNew()
	require.NotSame(t, &ct.value[0].Coeffs[0][0], &copyCT.value[0].Coeffs[0][0])
	copyCT.value[0].Coeffs[0][0]++
	require.Equal(t, before[0].Coeffs[0][0], ct.value[0].Coeffs[0][0])
	require.NoError(t, ct.ResizeDegree(2))
	require.Equal(t, 2, ct.Degree())
	require.Len(t, ct.value[2].Coeffs, 3)
	require.NoError(t, ct.ResizeDegree(0))
	require.Equal(t, 0, ct.Degree())
}

func TestFastStorageLevel0CoefficientImportWidths123(t *testing.T) {
	params := storageBoundaryParameters(t)
	source, expected := storageBoundarySource(t, params)
	for width := 1; width <= 3; width++ {
		fastCT, err := ImportLevel0(params, source, width, FastCiphertextDomain{})
		require.NoError(t, err, "storage width %d", width)
		require.Equal(t, 0, fastCT.LogicalLevel())
		require.Equal(t, width, fastCT.StorageWidth())
		require.Equal(t, source.Degree(), fastCT.Degree())
		require.True(t, source.Scale.Equal(fastCT.Scale()))
		requireFastLiftsEqual(t, fastCT, expected)
	}
}

func TestFastStorageLevel0NTTImportUsesCoefficientIntegers(t *testing.T) {
	params := storageBoundaryParameters(t)
	coeffSource, expected := storageBoundarySource(t, params)
	nttSource := coeffSource.CopyNew()
	logicalRing := params.RingQ().AtLevel(0)
	for component := range nttSource.Value {
		logicalRing.SubRings[0].NTT(nttSource.Value[component].Coeffs[0], nttSource.Value[component].Coeffs[0])
	}
	nttSource.IsNTT = true

	for width := 1; width <= 3; width++ {
		fastCT, err := ImportLevel0(params, nttSource, width, FastCiphertextDomain{IsNTT: true})
		require.NoError(t, err)
		require.True(t, fastCT.IsNTT())
		require.False(t, fastCT.IsMontgomery())
		requireFastLiftsEqual(t, fastCT, expected)
	}
}

func TestFastStorageExportLevel0AndHigherLogicalLevels(t *testing.T) {
	params := storageBoundaryParameters(t)
	source, expected := storageBoundarySource(t, params)
	fastCT, err := ImportLevel0(params, source, 3, FastCiphertextDomain{})
	require.NoError(t, err)

	level0, err := fastCT.ExportToLogical(FastCiphertextDomain{})
	require.NoError(t, err)
	require.Equal(t, 0, level0.Level())
	require.False(t, level0.IsNTT)
	require.False(t, level0.IsMontgomery)
	require.True(t, source.Scale.Equal(level0.Scale))
	for component := range source.Value {
		require.Equal(t, source.Value[component].Coeffs[0], level0.Value[component].Coeffs[0])
	}

	require.NoError(t, fastCT.SetLogicalLevel(2))
	higher, err := fastCT.ExportToLogical(FastCiphertextDomain{})
	require.NoError(t, err)
	require.Equal(t, 2, higher.Level())
	require.Equal(t, source.Degree(), higher.Degree())
	require.True(t, source.Scale.Equal(higher.Scale))
	fastParams := fastCT.Parameters()
	require.True(t, params.Equal(&fastParams))
	for component := range higher.Value {
		for k := 0; k < params.N(); k++ {
			x := expected[component][k]
			for logicalIndex := 0; logicalIndex <= 2; logicalIndex++ {
				q := params.Q()[logicalIndex]
				want := new(big.Int).Mod(new(big.Int).Set(x), new(big.Int).SetUint64(q)).Uint64()
				require.Equal(t, want, higher.Value[component].Coeffs[logicalIndex][k], "component=%d k=%d q-index=%d", component, k, logicalIndex)
			}
		}
	}
}

func TestFastStorageNTTToLogicalCoefficientAndNTTExport(t *testing.T) {
	params := storageBoundaryParameters(t)
	source, expected := storageBoundarySource(t, params)
	fastCT, err := ImportLevel0(params, source, 3, FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	require.NoError(t, fastCT.SetLogicalLevel(2))

	coeffOut, err := fastCT.ExportToLogical(FastCiphertextDomain{})
	require.NoError(t, err)
	require.False(t, coeffOut.IsNTT)
	nttOut, err := fastCT.ExportToLogical(FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	require.True(t, nttOut.IsNTT)
	require.False(t, nttOut.IsMontgomery)
	logicalRing := params.RingQ().AtLevel(2)
	for component := range coeffOut.Value {
		for logicalIndex := 0; logicalIndex <= 2; logicalIndex++ {
			wantNTT := make([]uint64, params.N())
			logicalRing.SubRings[logicalIndex].NTT(coeffOut.Value[component].Coeffs[logicalIndex], wantNTT)
			require.Equal(t, wantNTT, nttOut.Value[component].Coeffs[logicalIndex])
			for k := 0; k < params.N(); k++ {
				q := params.Q()[logicalIndex]
				want := new(big.Int).Mod(new(big.Int).Set(expected[component][k]), new(big.Int).SetUint64(q)).Uint64()
				require.Equal(t, want, coeffOut.Value[component].Coeffs[logicalIndex][k])
			}
		}
	}
}

func TestFastStorageConversionRejectsUnsupportedMontgomeryAndCapacity(t *testing.T) {
	params := storageBoundaryParameters(t)
	source, _ := storageBoundarySource(t, params)
	montgomerySource := source.CopyNew()
	montgomerySource.IsMontgomery = true
	_, err := ImportLevel0(params, montgomerySource, 1, FastCiphertextDomain{})
	require.ErrorContains(t, err, "Montgomery LogicalQ input")
	_, err = ImportLevel0(params, source, 1, FastCiphertextDomain{IsMontgomery: true})
	require.ErrorContains(t, err, "Montgomery Fast storage output")

	fastCT, err := ImportLevel0(params, source, 1, FastCiphertextDomain{})
	require.NoError(t, err)
	_, err = fastCT.ExportToLogical(FastCiphertextDomain{IsMontgomery: true})
	require.ErrorContains(t, err, "Montgomery LogicalQ output")
	fastCT.metadata.IsMontgomery = true
	_, err = fastCT.ExportToLogical(FastCiphertextDomain{})
	require.ErrorContains(t, err, "Montgomery Fast storage input")

	basis, err := newFastStorageBasis(params.LogN())
	require.NoError(t, err)
	capacity, err := basis.centeredCapacity(1)
	require.NoError(t, err)
	tooLarge := new(big.Int).Add(capacity, big.NewInt(1))
	_, err = encodeCenteredStorageValue(basis, storageInteger{magnitude: storageFromBig(tooLarge)}, 1)
	require.ErrorContains(t, err, "does not fit uniquely")

	for width := 1; width <= 3; width++ {
		actual, err := ImportLevel0(params, source, width, FastCiphertextDomain{})
		require.NoError(t, err, "actual logical q0 must fit width %d", width)
		require.Equal(t, width, actual.StorageWidth())
	}
}

func TestFastStorageLevel0ActualLogN13Q0CapacityWidths123(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 39},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	source, expected := storageBoundarySource(t, params)
	q0 := params.Q()[0]
	require.Equal(t, 55, bits.Len64(q0))

	basis, err := newFastStorageBasis(params.LogN())
	require.NoError(t, err)
	for width := 1; width <= 3; width++ {
		capacity, err := basis.centeredCapacity(width)
		require.NoError(t, err)
		require.True(t, new(big.Int).SetUint64(q0/2).Cmp(capacity) < 0, "actual q0 centered range must fit width %d", width)
		fastCT, err := ImportLevel0(params, source, width, FastCiphertextDomain{})
		require.NoError(t, err)
		requireFastLiftsEqual(t, fastCT, expected)
	}
}

func TestFastStorageConversionPreservesCKKSMetadata(t *testing.T) {
	params := storageBoundaryParameters(t)
	paramsBefore := params.ParametersLiteral()
	source, _ := storageBoundarySource(t, params)
	source.Scale = rlwe.NewScale(123456789)
	source.IsBatched = false
	source.LogDimensions = ring.Dimensions{Rows: 2, Cols: 3}
	source.IsBitReversed = true

	fastCT, err := ImportLevel0(params, source, 2, FastCiphertextDomain{})
	require.NoError(t, err)
	metadata := fastCT.MetaData()
	require.True(t, source.Scale.Equal(metadata.Scale))
	require.False(t, metadata.IsBatched)
	require.Equal(t, source.LogDimensions, metadata.LogDimensions)
	require.True(t, metadata.IsBitReversed)
	require.Equal(t, 0, fastCT.LogicalLevel())
	require.Equal(t, 2, fastCT.StorageWidth())
	require.True(t, params.Equal(&fastCT.params))

	output, err := fastCT.ExportToLogical(FastCiphertextDomain{})
	require.NoError(t, err)
	require.True(t, source.Scale.Equal(output.Scale))
	require.False(t, output.IsBatched)
	require.Equal(t, source.LogDimensions, output.LogDimensions)
	require.True(t, output.IsBitReversed)
	require.False(t, output.IsMontgomery)
	require.Equal(t, paramsBefore, params.ParametersLiteral(), "conversion must not mutate public CKKS parameters")
}

func storageBoundarySource(t *testing.T, params ckks.Parameters) (*rlwe.Ciphertext, [][]*big.Int) {
	t.Helper()
	source := ckks.NewCiphertext(params, 1, 0)
	source.Scale = rlwe.NewScale(1 << 30)
	source.IsNTT = false
	source.IsMontgomery = false
	source.IsBatched = true
	source.LogDimensions = params.LogMaxDimensions()
	q0 := params.Q()[0]
	values := []*big.Int{
		big.NewInt(0), big.NewInt(1), big.NewInt(-1), big.NewInt(12345), big.NewInt(-67890),
		new(big.Int).SetUint64(q0 / 2), new(big.Int).Neg(new(big.Int).SetUint64(q0 / 2)),
	}
	expected := make([][]*big.Int, len(source.Value))
	for component := range source.Value {
		expected[component] = make([]*big.Int, params.N())
		for k := 0; k < params.N(); k++ {
			value := new(big.Int).Set(values[(k+component)%len(values)])
			residue := new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(q0)).Uint64()
			source.Value[component].Coeffs[0][k] = residue
			centered := centeredStorageResidue(residue, q0)
			expected[component][k] = signedStorageBig(centered)
		}
	}
	return source, expected
}

func requireFastLiftsEqual(t *testing.T, fastCT *FastCiphertext, expected [][]*big.Int) {
	t.Helper()
	for component := range fastCT.value {
		coeffRows := make([][]uint64, fastCT.StorageWidth())
		for i := range coeffRows {
			coeffRows[i] = make([]uint64, fastCT.N())
			if fastCT.IsNTT() {
				subr, err := fastCT.basis.subring(i)
				require.NoError(t, err)
				subr.INTT(fastCT.value[component].Coeffs[i], coeffRows[i])
			} else {
				copy(coeffRows[i], fastCT.value[component].Coeffs[i])
			}
		}
		for k := 0; k < fastCT.N(); k++ {
			var residues [3]uint64
			for i := 0; i < fastCT.StorageWidth(); i++ {
				residues[i] = coeffRows[i][k]
			}
			value, err := fastCT.basis.decodeFixed(residues, fastCT.StorageWidth())
			require.NoError(t, err)
			require.Equal(t, expected[component][k], signedStorageBig(value), "component=%d coefficient=%d", component, k)
		}
	}
}

func cloneStorageRows(value []ring.Poly) []ring.Poly {
	copy := make([]ring.Poly, len(value))
	for i := range value {
		copy[i].Coeffs = make([][]uint64, len(value[i].Coeffs))
		for j := range value[i].Coeffs {
			copy[i].Coeffs[j] = append([]uint64(nil), value[i].Coeffs[j]...)
		}
	}
	return copy
}
