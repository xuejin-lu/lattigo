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

func q012TestParameters(t *testing.T) ckks.Parameters {
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
	require.Equal(t, 56, bits.Len64(params.Q()[0]))
	require.LessOrEqual(t, bits.Len64(params.Q()[1]), 39)
	require.LessOrEqual(t, bits.Len64(params.Q()[2]), 40)
	return params
}

func uint192Big(value uint192) *big.Int {
	result := new(big.Int).SetUint64(value.hi)
	result.Lsh(result, 64)
	result.Add(result, new(big.Int).SetUint64(value.mid))
	result.Lsh(result, 64)
	result.Add(result, new(big.Int).SetUint64(value.lo))
	return result
}

func bigUint192(value *big.Int) uint192 {
	mask := new(big.Int).Lsh(big.NewInt(1), 64)
	mask.Sub(mask, big.NewInt(1))
	lo := new(big.Int).And(new(big.Int).Set(value), mask).Uint64()
	value = new(big.Int).Rsh(new(big.Int).Set(value), 64)
	mid := new(big.Int).And(new(big.Int).Set(value), mask).Uint64()
	hi := new(big.Int).Rsh(value, 64).Uint64()
	return uint192{lo: lo, mid: mid, hi: hi}
}

func TestQ012FixedWidthCRTAndCenteredReconstruction(t *testing.T) {
	params := q012TestParameters(t)
	q := params.Q()
	q0, q1, q2 := q[0], q[1], q[2]
	q0InvQ1, ok := inverseMod(q0%q1, q1)
	require.True(t, ok)
	q01Lo, q01Hi := bitsMul64ForTest(q0, q1)
	q01InvQ2, ok := inverseMod(mod128By64(q01Lo, q01Hi, q2), q2)
	require.True(t, ok)
	modulus, half, _, _ := q012Modulus(q0, q1, q2)
	Q := new(big.Int).Mul(new(big.Int).SetUint64(q0), new(big.Int).SetUint64(q1))
	Q.Mul(Q, new(big.Int).SetUint64(q2))

	values := []*big.Int{
		big.NewInt(0), big.NewInt(1), big.NewInt(-1),
		new(big.Int).Lsh(big.NewInt(1), 80),
		new(big.Int).Neg(new(big.Int).Lsh(big.NewInt(1), 80)),
		new(big.Int).Sub(new(big.Int).Rsh(new(big.Int).Set(Q), 1), big.NewInt(1)),
		new(big.Int).Add(new(big.Int).Rsh(new(big.Int).Set(Q), 1), big.NewInt(1)),
	}
	for _, signed := range values {
		residue := func(modulus uint64) uint64 {
			return new(big.Int).Mod(new(big.Int).Set(signed), new(big.Int).SetUint64(modulus)).Uint64()
		}
		got := crtQ012(residue(q0), residue(q1), residue(q2), q0, q1, q2, q0InvQ1, q01InvQ2)
		wantResidue := new(big.Int).Mod(new(big.Int).Set(signed), Q)
		require.Equal(t, wantResidue, uint192Big(got), "CRT value %s", signed)
		wantNegative := wantResidue.Cmp(new(big.Int).Rsh(new(big.Int).Set(Q), 1)) > 0
		wantMagnitude := new(big.Int).Set(wantResidue)
		if wantNegative {
			wantMagnitude.Sub(Q, wantMagnitude)
		}
		gotMagnitude, gotNegative := centeredQ012(got, modulus, half)
		require.Equal(t, wantNegative, gotNegative, "center sign %s", signed)
		require.Equal(t, wantMagnitude, uint192Big(gotMagnitude), "center magnitude %s", signed)
	}
}

func TestQ012FixedWidthRoundedDivideMatchesBigInt(t *testing.T) {
	params := q012TestParameters(t)
	q := params.Q()
	Q := new(big.Int).Mul(new(big.Int).SetUint64(q[0]), new(big.Int).SetUint64(q[1]))
	Q.Mul(Q, new(big.Int).SetUint64(q[2]))
	values := []*big.Int{
		big.NewInt(0), big.NewInt(1), big.NewInt(2),
		new(big.Int).Sub(new(big.Int).Rsh(new(big.Int).Set(Q), 1), big.NewInt(1)),
		new(big.Int).Lsh(big.NewInt(1), 100),
	}
	for _, value := range values {
		for _, divisor := range []uint64{q[3], q[2], uint64(1) << 32} {
			got := roundedMagnitude192(bigUint192(value), divisor)
			want, remainder := new(big.Int).QuoRem(new(big.Int).Set(value), new(big.Int).SetUint64(divisor), new(big.Int))
			if new(big.Int).Lsh(remainder, 1).Cmp(new(big.Int).SetUint64(divisor)) > 0 {
				want.Add(want, big.NewInt(1))
			}
			require.True(t, want.Cmp(uint192Big(got)) == 0, "value=%s divisor=%d want=%s got=%s", value, divisor, want, uint192Big(got))
		}
	}
}

func TestQ012FastRescaleMatchesStandardForBoundedCoefficients(t *testing.T) {
	params := q012TestParameters(t)
	fastEval := NewEvaluator(params)
	standardEval := ckks.NewEvaluator(params, nil)
	fast := NewCiphertext(params, 1, 3)
	standard := ckks.NewCiphertext(params, 1, 3)
	fast.IsNTT, standard.IsNTT = true, true
	fast.Scale, standard.Scale = rlwe.NewScale(1<<30), rlwe.NewScale(1<<30)

	values := []int64{0, 1, -1, 123, -456, 1 << 20, -(1 << 20)}
	for component := 0; component <= 1; component++ {
		for i := range values {
			x := big.NewInt(values[i])
			for limb, modulus := range params.Q() {
				residue := new(big.Int).Mod(new(big.Int).Set(x), new(big.Int).SetUint64(modulus)).Uint64()
				if limb < 3 {
					fast.Value[component].Coeffs[limb][i] = residue
				}
				standard.Value[component].Coeffs[limb][i] = residue
			}
		}
		for limb := 0; limb <= 2; limb++ {
			params.RingQ().SubRings[limb].NTT(fast.Value[component].Coeffs[limb], fast.Value[component].Coeffs[limb])
		}
		for limb := 0; limb <= 3; limb++ {
			params.RingQ().SubRings[limb].NTT(standard.Value[component].Coeffs[limb], standard.Value[component].Coeffs[limb])
		}
	}
	fastOut := NewCiphertext(params, 1, 2)
	standardOut := ckks.NewCiphertext(params, 1, 2)
	fastOut.IsNTT, standardOut.IsNTT = true, true
	require.NoError(t, fastEval.Rescale(fast, fastOut))
	require.NoError(t, standardEval.Rescale(standard, standardOut))
	for component := 0; component <= 1; component++ {
		for limb := 0; limb <= 2; limb++ {
			require.Equal(t, standardOut.Value[component].Coeffs[limb], fastOut.Value[component].Coeffs[limb], "component=%d limb=%d", component, limb)
		}
	}
}

func bitsMul64ForTest(a, b uint64) (uint64, uint64) {
	hi, lo := bits.Mul64(a, b)
	return lo, hi
}
