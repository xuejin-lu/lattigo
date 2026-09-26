package fast

import (
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFastCompactCiphertextStorage(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 39, 39, 39, 39, 39, 39, 39},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)

	for _, level := range []int{0, 1, 2, 3, 4, 8} {
		ct := NewCiphertext(params, 1, level)
		requireQPrefixStorage(t, ct, level, params.N())
	}
}

func TestFastCiphertextResizeIsStructuralAndPreservesMetadata(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{50, 50, 50, 50, 50, 50, 50, 50, 50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)

	ct := NewCiphertext(params, 1, 3)
	for component := range ct.Value {
		for limb := 0; limb < 4; limb++ {
			for i := range ct.Value[component].Coeffs[limb] {
				ct.Value[component].Coeffs[limb][i] = uint64(100*component + 10*limb + i + 1)
			}
		}
	}
	ct.IsNTT = true
	ct.IsMontgomery = true
	ct.Scale = rlwe.NewScale(1 << 30)
	metadata := ct.MetaData.CopyNew()
	originalRows := make([][][]uint64, len(ct.Value))
	for component := range ct.Value {
		originalRows[component] = make([][]uint64, 4)
		for limb := 0; limb < 4; limb++ {
			originalRows[component][limb] = append([]uint64(nil), ct.Value[component].Coeffs[limb]...)
		}
	}

	for _, level := range []int{2, 1, 0} {
		Resize(ct, 1, level, params.N())
		requireQPrefixStorage(t, ct, level, params.N())
		for component := range ct.Value {
			for limb := 0; limb < level+1; limb++ {
				require.Equal(t, originalRows[component][limb], ct.Value[component].Coeffs[limb])
			}
		}
	}
	require.Equal(t, metadata, ct.MetaData)
	require.True(t, ct.IsNTT)
	require.True(t, ct.IsMontgomery)
	require.Equal(t, rlwe.NewScale(1<<30), ct.Scale)

	Resize(ct, 1, 3, params.N())
	requireQPrefixStorage(t, ct, 3, params.N())
	for component := range ct.Value {
		require.Equal(t, originalRows[component][0], ct.Value[component].Coeffs[0])
		for limb := 1; limb < 4; limb++ {
			for _, coefficient := range ct.Value[component].Coeffs[limb] {
				require.Zero(t, coefficient, "grown row %d must not infer or restore a representative", limb)
			}
		}
	}
}

func TestFastCiphertextCopyAndDegreeGrowthPreservePrefixShape(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{50, 50, 50, 50, 50, 50, 50, 50, 50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)

	src := NewCiphertext(params, 1, 8)
	for component := range src.Value {
		for limb := 0; limb < 4; limb++ {
			for i := range src.Value[component].Coeffs[limb] {
				src.Value[component].Coeffs[limb][i] = uint64(1000*component + 100*limb + i + 1)
			}
		}
	}
	src.IsNTT = true
	src.IsMontgomery = true
	src.Scale = rlwe.NewScale(1 << 30)

	copyDst := NewCiphertext(params, 1, 8)
	copyDst.Copy(src)
	copyNew := src.CopyNew()
	for _, copied := range []*rlwe.Ciphertext{copyDst, copyNew} {
		requireQPrefixStorage(t, copied, 8, params.N())
		require.Equal(t, src.MetaData, copied.MetaData)
		for component := range src.Value {
			for limb := 0; limb < 4; limb++ {
				require.Equal(t, src.Value[component].Coeffs[limb], copied.Value[component].Coeffs[limb])
				require.NotSame(t, &src.Value[component].Coeffs[limb][0], &copied.Value[component].Coeffs[limb][0])
			}
		}
	}

	Resize(src, 3, 8, params.N())
	requireQPrefixStorage(t, src, 8, params.N())
	require.Equal(t, 3, src.Degree())
	for component := 2; component <= 3; component++ {
		for limb := 0; limb < 4; limb++ {
			require.Len(t, src.Value[component].Coeffs[limb], params.N())
			for _, coefficient := range src.Value[component].Coeffs[limb] {
				require.Zero(t, coefficient)
			}
		}
	}
	for component := 0; component <= 1; component++ {
		for limb := 0; limb < 4; limb++ {
			require.Equal(t, copyNew.Value[component].Coeffs[limb], src.Value[component].Coeffs[limb])
		}
	}
}

func TestFastResizeCompactsFullCiphertextBacking(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{50, 50, 50, 50, 50, 50, 50, 50, 50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	ct := ckks.NewCiphertext(params, 1, 8)
	for limb := range ct.Value[0].Coeffs {
		for i := range ct.Value[0].Coeffs[limb] {
			ct.Value[0].Coeffs[limb][i] = uint64(limb + i + 1)
		}
	}
	retained := make([][]uint64, 4)
	for limb := range retained {
		retained[limb] = append([]uint64(nil), ct.Value[0].Coeffs[limb]...)
	}

	Resize(ct, 1, 8, params.N())
	requireQPrefixStorage(t, ct, 8, params.N())
	for limb := range retained {
		require.Equal(t, retained[limb], ct.Value[0].Coeffs[limb])
	}
}

func TestFastEvaluatorScratchUsesQPrefixWidth(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	width, err := QPrefixWidth(params.MaxLevel())
	require.NoError(t, err)

	for _, poly := range eval.nttScratch {
		require.Len(t, poly.Coeffs, width)
		for limb := 0; limb < width; limb++ {
			require.Len(t, poly.Coeffs[limb], params.N())
		}
	}
	for _, poly := range [][][]uint64{
		eval.linearTransformScratch.acc0.Coeffs,
		eval.linearTransformScratch.acc1.Coeffs,
		eval.linearTransformScratch.inner0.Coeffs,
		eval.linearTransformScratch.inner1.Coeffs,
		eval.linearTransformScratch.outer0.Coeffs,
		eval.linearTransformScratch.outer1.Coeffs,
		eval.linearTransformScratch.term0.Coeffs,
		eval.linearTransformScratch.term1.Coeffs,
		eval.linearTransformScratch.rot0.Coeffs,
		eval.linearTransformScratch.rot1.Coeffs,
	} {
		require.Len(t, poly, width)
	}
	require.Len(t, eval.automorphismScratch, width)
	for _, scratch := range eval.automorphismScratch {
		require.Len(t, scratch, params.N())
	}
	require.Len(t, eval.rescaleScratch.coeff.Coeffs, width)
	require.Len(t, eval.rescaleScratch.result.Coeffs, width)
}

func TestFastAliasedAddPreservesLegacyUnwrittenPrefixRows(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	ct := NewCiphertext(params, 1, 3)
	ct.Scale = params.DefaultScale()
	ct.IsNTT = true
	for component := range ct.Value {
		for i := range ct.Value[component].Coeffs[0] {
			ct.Value[component].Coeffs[0][i] = 7
			ct.Value[component].Coeffs[1][i] = 11
			ct.Value[component].Coeffs[2][i] = 19
			ct.Value[component].Coeffs[3][i] = 23
		}
	}

	require.NoError(t, eval.Add(ct, ct, ct))
	requireQPrefixStorage(t, ct, 3, params.N())
	for component := range ct.Value {
		require.Equal(t, uint64(14), ct.Value[component].Coeffs[0][0])
		require.Equal(t, uint64(22), ct.Value[component].Coeffs[1][0])
		require.Equal(t, uint64(19), ct.Value[component].Coeffs[2][0])
		require.Equal(t, uint64(23), ct.Value[component].Coeffs[3][0])
	}
}

func TestFastPointwiseScratchDoesNotPromoteUnwrittenPrefixRows(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	a := NewCiphertext(params, 1, 3)
	b := NewCiphertext(params, 1, 3)
	out := NewCiphertext(params, 1, 3)
	for component := range a.Value {
		for i := range a.Value[component].Coeffs[0] {
			for limb := 0; limb < 2; limb++ {
				a.Value[component].Coeffs[limb][i] = 2
				b.Value[component].Coeffs[limb][i] = 3
			}
			out.Value[component].Coeffs[2][i] = 19
			out.Value[component].Coeffs[3][i] = 23
		}
	}

	require.NoError(t, eval.pointMul(a.Value[0], b.Value[0], out.Value[0], false))
	for limb := 0; limb < 2; limb++ {
		q := params.Q()[limb]
		for _, coefficient := range out.Value[0].Coeffs[limb] {
			require.Equal(t, uint64(6)%q, coefficient)
		}
	}
	require.Equal(t, uint64(19), out.Value[0].Coeffs[2][0])
	require.Equal(t, uint64(23), out.Value[0].Coeffs[3][0])

	for i := range out.Value[0].Coeffs[0] {
		out.Value[0].Coeffs[0][i] = 5
		out.Value[0].Coeffs[1][i] = 7
	}
	require.NoError(t, eval.pointMulThenAdd(a.Value[0], b.Value[0], out.Value[0], false))
	require.Equal(t, uint64(11), out.Value[0].Coeffs[0][0])
	require.Equal(t, uint64(13), out.Value[0].Coeffs[1][0])
	require.Equal(t, uint64(19), out.Value[0].Coeffs[2][0])
	require.Equal(t, uint64(23), out.Value[0].Coeffs[3][0])
}

func requireQPrefixStorage(t *testing.T, ct *rlwe.Ciphertext, level, n int) {
	t.Helper()
	width, err := QPrefixWidth(level)
	require.NoError(t, err)
	require.Equal(t, level, ct.Level())
	require.Equal(t, n, ct.N())
	for component := range ct.Value {
		require.Len(t, ct.Value[component].Coeffs, level+1)
		for limb, row := range ct.Value[component].Coeffs {
			if limb < width {
				require.Len(t, row, n, "component=%d limb=%d", component, limb)
			} else {
				require.Empty(t, row, "dormant component=%d limb=%d must have no backing", component, limb)
				require.Zero(t, cap(row), "dormant component=%d limb=%d must have zero capacity", component, limb)
			}
		}
	}
}

func TestFastCompactCiphertextArithmeticAndPublicLevelOneStorage(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 39, 39},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewEvaluator(params)
	left := NewCiphertext(params, 1, 3)
	right := NewCiphertext(params, 1, 3)
	out := NewCiphertext(params, 1, 3)
	for _, ct := range []*rlwe.Ciphertext{left, right, out} {
		ct.IsNTT = true
		ct.IsMontgomery = true
		ct.Scale = params.DefaultScale()
	}
	for d := range left.Value {
		for limb := 0; limb < 2; limb++ {
			for i := range left.Value[d].Coeffs[limb] {
				left.Value[d].Coeffs[limb][i] = uint64(i+1+d+limb) % params.RingQ().SubRings[limb].Modulus
				right.Value[d].Coeffs[limb][i] = uint64(2*i+3+d+limb) % params.RingQ().SubRings[limb].Modulus
			}
		}
	}
	require.NoError(t, eval.Add(left, right, out))
	for d := range out.Value {
		requireQPrefixStorage(t, out, out.Level(), params.N())
		for limb := maintainedLimbCount(&params, out.Level()); limb <= out.Level(); limb++ {
			for _, coefficient := range out.Value[d].Coeffs[limb] {
				require.Zero(t, coefficient, "unwritten prefix limb %d must not be treated as authoritative", limb)
			}
		}
	}

	levelOne := NewCiphertext(params, 1, 1)
	require.Len(t, levelOne.Value[0].Coeffs[0], params.N())
	require.Len(t, levelOne.Value[0].Coeffs[1], params.N())
}
