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
		LogQ:            []int{55, 39, 39, 39, 39},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)

	ct := NewCiphertext(params, 1, 4)
	require.Equal(t, 4, ct.Level())
	require.Equal(t, params.N(), ct.N())
	require.Equal(t, 1, ct.Degree())
	for d := range ct.Value {
		require.Len(t, ct.Value[d].Coeffs, 5)
		require.Len(t, ct.Value[d].Coeffs[0], params.N())
		require.Len(t, ct.Value[d].Coeffs[1], params.N())
		for limb := 2; limb <= ct.Level(); limb++ {
			require.Empty(t, ct.Value[d].Coeffs[limb])
		}
	}

	Resize(ct, 1, 0, params.N())
	require.Equal(t, 0, ct.Level())
	for d := range ct.Value {
		require.Len(t, ct.Value[d].Coeffs, 1)
		require.Len(t, ct.Value[d].Coeffs[0], params.N())
	}

	Resize(ct, 1, 4, params.N())
	require.Equal(t, 4, ct.Level())
	for d := range ct.Value {
		require.Len(t, ct.Value[d].Coeffs[0], params.N())
		require.Len(t, ct.Value[d].Coeffs[1], params.N())
		for limb := 2; limb <= ct.Level(); limb++ {
			require.Empty(t, ct.Value[d].Coeffs[limb])
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
		for limb := 0; limb < 2; limb++ {
			require.Len(t, out.Value[d].Coeffs[limb], params.N())
		}
		for limb := 2; limb <= out.Level(); limb++ {
			require.Empty(t, out.Value[d].Coeffs[limb])
		}
	}

	levelOne := NewCiphertext(params, 1, 1)
	require.Len(t, levelOne.Value[0].Coeffs[0], params.N())
	require.Len(t, levelOne.Value[0].Coeffs[1], params.N())
}
