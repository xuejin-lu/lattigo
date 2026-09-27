package fast

import (
	"fmt"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestDropLevelSameLiftFitsStrictTargetPrefix(t *testing.T) {
	params := rescaleTestParameters(t)
	eval := NewEvaluator(params)
	for _, transition := range []struct{ source, target int }{{3, 2}, {2, 1}, {1, 0}} {
		t.Run(fmt.Sprintf("%d_to_%d", transition.source, transition.target), func(t *testing.T) {
			rows, err := QPrefixWidth(transition.target)
			require.NoError(t, err)
			targetProduct, err := QPrefixProduct(params.Q(), transition.target)
			require.NoError(t, err)
			half := new(big.Int).Rsh(targetProduct, 1)
			for _, representation := range []struct {
				name         string
				isNTT        bool
				isMontgomery bool
			}{{"coeff", false, false}, {"ntt", true, false}, {"ntt_montgomery", true, true}} {
				t.Run(representation.name, func(t *testing.T) {
					for _, sign := range []int64{1, -1} {
						value := new(big.Int).Mul(new(big.Int).Set(half), big.NewInt(sign))
						in := makeLevelTransitionCiphertext(params, transition.source, value, representation.isNTT, representation.isMontgomery)
						out := NewCiphertext(params, in.Degree(), transition.target)
						require.NoError(t, eval.DropLevelSameLift(in, transition.target, out))
						require.Equal(t, transition.target, out.Level())
						require.Equal(t, in.IsNTT, out.IsNTT)
						require.Equal(t, in.IsMontgomery, out.IsMontgomery)
						require.Equal(t, *in.MetaData, *out.MetaData)
						require.True(t, in.Scale.Equal(out.Scale))
						for component := range out.Value {
							for row := 0; row < rows; row++ {
								require.Equal(t, in.Value[component].Coeffs[row], out.Value[component].Coeffs[row])
							}
						}
					}
				})
			}
		})
	}
}

func TestDropLevelSameLiftRejectsFirstNonFittingValueTransactionally(t *testing.T) {
	params := rescaleTestParameters(t)
	eval := NewEvaluator(params)
	target := 2
	targetProduct, err := QPrefixProduct(params.Q(), target)
	require.NoError(t, err)
	firstOutside := new(big.Int).Add(new(big.Int).Rsh(targetProduct, 1), big.NewInt(1))
	in := makeLevelTransitionCiphertext(params, 3, firstOutside, true, true)
	out := NewCiphertext(params, 1, target)
	fillFastCKKSCiphertext(out, params, [][]int64{{17, -19}, {23, -29}}, 41)
	out.IsNTT, out.IsMontgomery = false, false
	out.Scale = rlwe.NewScale(12345)
	inBefore, outBefore := in.CopyNew(), out.CopyNew()
	err = eval.DropLevelSameLift(in, target, out)
	var capacityErr *QPrefixCapacityError
	require.ErrorAs(t, err, &capacityErr)
	require.Equal(t, target, capacityErr.Level)
	require.Equal(t, inBefore.Level(), in.Level())
	require.Equal(t, *inBefore.MetaData, *in.MetaData)
	require.Equal(t, inBefore.Value, in.Value)
	require.Equal(t, *outBefore.MetaData, *out.MetaData)
	require.Equal(t, outBefore.Value, out.Value)

	inPlaceBefore := in.CopyNew()
	err = eval.DropLevelSameLift(in, target, in)
	require.ErrorAs(t, err, &capacityErr)
	require.Equal(t, *inPlaceBefore.MetaData, *in.MetaData)
	require.Equal(t, inPlaceBefore.Value, in.Value)
}

func TestDropLevelCanonicalUsesTargetCenteredRepresentative(t *testing.T) {
	params := rescaleTestParameters(t)
	eval := NewEvaluator(params)
	target := 2
	targetProduct, err := QPrefixProduct(params.Q(), target)
	require.NoError(t, err)
	half := new(big.Int).Rsh(new(big.Int).Set(targetProduct), 1)
	sourceLift := new(big.Int).Add(new(big.Int).Set(half), big.NewInt(1))
	canonical := new(big.Int).Neg(new(big.Int).Set(half))
	require.NotEqual(t, sourceLift, canonical)
	in := makeLevelTransitionCiphertext(params, 3, sourceLift, false, false)
	out := NewCiphertext(params, 1, target)
	require.NoError(t, eval.DropLevelCanonical(in, target, out))
	require.Equal(t, target, out.Level())
	require.True(t, in.Scale.Equal(out.Scale))
	for component := range out.Value {
		for row := 0; row <= target; row++ {
			wantResidue := new(big.Int).Mod(new(big.Int).Set(canonical), new(big.Int).SetUint64(params.Q()[row])).Uint64()
			for coefficient, got := range out.Value[component].Coeffs[row] {
				require.Equal(t, wantResidue, got, "component=%d q%d coefficient=%d", component, row, coefficient)
			}
		}
	}
	var capacityErr *QPrefixCapacityError
	require.ErrorAs(t, eval.DropLevelSameLift(in, target, NewCiphertext(params, 1, target)), &capacityErr)
}

func makeLevelTransitionCiphertext(params ckks.Parameters, level int, value *big.Int, isNTT, isMontgomery bool) *rlwe.Ciphertext {
	ct := NewCiphertext(params, 1, level)
	ct.IsNTT = isNTT
	ct.IsMontgomery = isMontgomery
	ct.Scale = rlwe.NewScale(1 << 40)
	rows, _ := QPrefixWidth(level)
	for component := range ct.Value {
		for row := 0; row < rows; row++ {
			q := params.Q()[row]
			residue := new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(q)).Uint64()
			for coefficient := range ct.Value[component].Coeffs[row] {
				ct.Value[component].Coeffs[row][coefficient] = residue
			}
			if isNTT {
				params.RingQ().SubRings[row].NTT(ct.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
			}
			if isMontgomery {
				params.RingQ().SubRings[row].MForm(ct.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
			}
		}
	}
	return ct
}
