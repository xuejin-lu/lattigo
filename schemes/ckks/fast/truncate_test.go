package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFastTruncateDegree2To1(t *testing.T) {
	params := testFastCKKSParameters(t)
	for _, level := range []int{1, 2, 3} {
		in := ckks.NewCiphertext(params, 2, level)
		out := ckks.NewCiphertext(params, 2, level)
		in.IsNTT = true
		in.IsMontgomery = true
		in.IsBatched = true
		in.Scale = rlwe.NewScale(big.NewFloat(123.5))
		in.LogDimensions.Rows = 2
		in.LogDimensions.Cols = 3
		fillFastCKKSCiphertext(in, params, nil, uint64(level+5))
		fillFastCKKSCiphertext(out, params, nil, uint64(level+50))

		want0 := append([][]uint64(nil), in.Value[0].Coeffs[0:2]...)
		want1 := append([][]uint64(nil), in.Value[1].Coeffs[0:2]...)
		dormant := cloneDormant(out, 2)
		metadata := *in.MetaData

		require.NoError(t, FastTruncateDegree2To1(in, out))
		require.Equal(t, 1, out.Degree())
		require.Equal(t, level, out.Level())
		require.Equal(t, metadata, *out.MetaData)
		require.Equal(t, want0[0], out.Value[0].Coeffs[0])
		require.Equal(t, want0[1], out.Value[0].Coeffs[1])
		require.Equal(t, want1[0], out.Value[1].Coeffs[0])
		require.Equal(t, want1[1], out.Value[1].Coeffs[1])
		// c0/c1 dormant limbs are not written by truncation.
		require.Equal(t, dormant[:2], cloneDormant(out, 2))
	}
}

func TestFastTruncateDegree2To1IgnoresC2(t *testing.T) {
	params := testFastCKKSParameters(t)
	base := ckks.NewCiphertext(params, 2, 3)
	fillFastCKKSCiphertext(base, params, nil, 7)

	outs := make([]*rlwe.Ciphertext, 3)
	for i := range outs {
		in := base.CopyNew()
		for j := range in.Value[2].Coeffs {
			for k := range in.Value[2].Coeffs[j] {
				in.Value[2].Coeffs[j][k] = uint64(i+1)*0x9e3779b9 + uint64(j+k)
			}
		}
		outs[i] = ckks.NewCiphertext(params, 1, 3)
		require.NoError(t, FastTruncateDegree2To1(in, outs[i]))
	}
	for i := 1; i < len(outs); i++ {
		require.Equal(t, outs[0].Value[0].Coeffs[:2], outs[i].Value[0].Coeffs[:2])
		require.Equal(t, outs[0].Value[1].Coeffs[:2], outs[i].Value[1].Coeffs[:2])
	}
}

func TestFastTruncateDegree2To1InPlace(t *testing.T) {
	params := testFastCKKSParameters(t)
	ct := ckks.NewCiphertext(params, 2, 3)
	fillFastCKKSCiphertext(ct, params, nil, 91)
	want0 := append([][]uint64(nil), ct.Value[0].Coeffs[0:2]...)
	want1 := append([][]uint64(nil), ct.Value[1].Coeffs[0:2]...)

	require.NoError(t, FastTruncateDegree2To1(ct, ct))
	require.Equal(t, 1, ct.Degree())
	require.Equal(t, want0[0], ct.Value[0].Coeffs[0])
	require.Equal(t, want0[1], ct.Value[0].Coeffs[1])
	require.Equal(t, want1[0], ct.Value[1].Coeffs[0])
	require.Equal(t, want1[1], ct.Value[1].Coeffs[1])
}

func TestFastTruncateDegree2To1RejectsInvalidInputs(t *testing.T) {
	params := testFastCKKSParameters(t)
	degree1 := ckks.NewCiphertext(params, 1, 3)
	degree2 := ckks.NewCiphertext(params, 2, 3)
	wrongLevel := ckks.NewCiphertext(params, 1, 2)
	for _, tc := range []struct {
		name string
		in   *rlwe.Ciphertext
		out  *rlwe.Ciphertext
	}{
		{"nil input", nil, degree1},
		{"nil output", degree2, nil},
		{"input degree", degree1, degree1},
		{"output degree", degree2, ckks.NewCiphertext(params, 0, 3)},
		{"level mismatch", degree2, wrongLevel},
	} {
		t.Run(tc.name, func(t *testing.T) {
			require.Error(t, FastTruncateDegree2To1(tc.in, tc.out))
		})
	}
}
