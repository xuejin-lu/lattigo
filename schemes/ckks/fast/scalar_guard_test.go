package fast

import (
	"math/big"
	"math/bits"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func TestRoundedMagnitude128ByTwoMatchesBigInt(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	q0, q1 := params.Q()[0], params.Q()[1]
	q := new(big.Int).Mul(new(big.Int).SetUint64(q0), new(big.Int).SetUint64(q1))
	half := new(big.Int).Rsh(new(big.Int).Set(q), 1)
	values := []*big.Int{
		new(big.Int), big.NewInt(1), big.NewInt(2), big.NewInt(3),
		big.NewInt(-1), big.NewInt(-2), big.NewInt(-3),
		new(big.Int).Sub(new(big.Int).Set(half), big.NewInt(1)),
		new(big.Int).Neg(new(big.Int).Sub(new(big.Int).Set(half), big.NewInt(1))),
	}
	qHi, qLo := bits.Mul64(q0, q1)
	halfLo := (qLo >> 1) | (qHi << 63)
	halfHi := qHi >> 1
	inverse, ok := inverseMod(q0%q1, q1)
	require.True(t, ok)
	for _, value := range values {
		r0 := new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(q0)).Uint64()
		r1 := new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(q1)).Uint64()
		xLo, xHi := crtQ01(r0, r1, q0, q1, inverse)
		magnitudeLo, magnitudeHi, negative := roundedMagnitude128ByTwo(xLo, xHi, qLo, qHi, halfLo, halfHi)
		got := new(big.Int).SetBits([]big.Word{big.Word(magnitudeLo), big.Word(magnitudeHi)})
		if negative {
			got.Neg(got)
		}
		want := new(big.Int).Set(value)
		negativeWant := want.Sign() < 0
		if negativeWant {
			want.Abs(want)
		}
		want.Add(want, big.NewInt(1)).Rsh(want, 1)
		if negativeWant {
			want.Neg(want)
		}
		require.Equal(t, want.String(), got.String(), "value=%s", value)
	}
}

func TestCenteredRoundedDivideByTwoPreservesBothComponentsAndMontgomery(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewEvaluator(params)
	ct := NewCiphertext(params, 1, 1)
	ct.IsNTT, ct.IsMontgomery = true, true
	ct.Scale = rlwe.NewScale(1 << 30)
	q0, q1 := params.Q()[0], params.Q()[1]
	values := [][]int64{
		{0, 1, -1, 2, -2, 3, -3, 5, -5, 7, -7, 0, 1, -1, 2, -2},
		{11, -11, 13, -13, 17, -17, 19, -19, 23, -23, 29, -29, 31, -31, 37, -37},
	}
	for component, row := range values {
		coeff := params.RingQ().AtLevel(1).NewPoly()
		for i, value := range row {
			for limb, modulus := range []uint64{q0, q1} {
				coeff.Coeffs[limb][i] = new(big.Int).Mod(new(big.Int).SetInt64(value), new(big.Int).SetUint64(modulus)).Uint64()
			}
		}
		for limb := 0; limb < 2; limb++ {
			params.RingQ().SubRings[limb].NTT(coeff.Coeffs[limb], ct.Value[component].Coeffs[limb])
			params.RingQ().SubRings[limb].MForm(ct.Value[component].Coeffs[limb], ct.Value[component].Coeffs[limb])
		}
	}
	original := *ct.MetaData
	require.NoError(t, eval.contractCenteredRoundedDivideByTwo(ct))
	require.Equal(t, original, *ct.MetaData)
	require.True(t, ct.IsNTT)
	require.True(t, ct.IsMontgomery)
	for component, row := range values {
		got := params.RingQ().AtLevel(1).NewPoly()
		for limb := 0; limb < 2; limb++ {
			params.RingQ().SubRings[limb].IMForm(ct.Value[component].Coeffs[limb], got.Coeffs[limb])
			params.RingQ().SubRings[limb].INTT(got.Coeffs[limb], got.Coeffs[limb])
		}
		for i, value := range row {
			want := new(big.Int).SetInt64(value)
			negative := want.Sign() < 0
			if negative {
				want.Abs(want)
			}
			want.Add(want, big.NewInt(1)).Rsh(want, 1)
			if negative {
				want.Neg(want)
			}
			for limb, modulus := range []uint64{q0, q1} {
				wantResidue := new(big.Int).Mod(new(big.Int).Set(want), new(big.Int).SetUint64(modulus)).Uint64()
				require.Equal(t, wantResidue, got.Coeffs[limb][i], "component=%d coeff=%d limb=%d", component, i, limb)
			}
		}
	}
}

func TestMulThenAddOneBitScalarGuardRestoresNativeMetadata(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewEvaluator(params)
	power := NewCiphertext(params, 1, 1)
	power.IsNTT, power.IsMontgomery = true, true
	power.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	accumulator := NewCiphertext(params, 1, 1)
	accumulator.IsNTT, accumulator.IsMontgomery = true, true
	accumulator.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 92))
	native := accumulator.CopyNew()
	coefficient := bignum.ToComplex(0.5, 128)
	require.NoError(t, eval.MulThenAdd(power, coefficient, native))
	nativeScale := native.Scale
	nativeMeta := *native.MetaData
	guarded := accumulator.CopyNew()
	require.NoError(t, eval.MulThenAddOneBitScalarGuard(power, coefficient, guarded))
	require.Equal(t, native.Level(), guarded.Level())
	require.Equal(t, native.Degree(), guarded.Degree())
	require.Equal(t, native.IsNTT, guarded.IsNTT)
	require.Equal(t, native.IsMontgomery, guarded.IsMontgomery)
	require.True(t, nativeScale.Equal(guarded.Scale))
	require.Equal(t, nativeMeta, *guarded.MetaData)
	for d := range native.Value {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, native.Value[d].Coeffs[limb], guarded.Value[d].Coeffs[limb])
		}
	}
}
