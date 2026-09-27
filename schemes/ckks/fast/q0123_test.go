package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestQ0123FixedWidthProductCRTAndCenterMatchBigInt(t *testing.T) {
	params := q012TestParameters(t)
	q := params.Q()
	q0, q1, q2, q3 := q[0], q[1], q[2], q[3]
	modulus, half, q012, overflow := q0123Modulus(q0, q1, q2, q3)
	require.False(t, overflow)
	wantModulus := new(big.Int).SetUint64(q0)
	for _, qi := range q[1:4] {
		wantModulus.Mul(wantModulus, new(big.Int).SetUint64(qi))
	}
	wantHalf := new(big.Int).Rsh(new(big.Int).Set(wantModulus), 1)
	require.Equal(t, wantModulus, uint192Big(modulus))
	require.Equal(t, wantHalf, uint192Big(half))
	wantQ012 := new(big.Int).Quo(new(big.Int).Set(wantModulus), new(big.Int).SetUint64(q3))
	require.Equal(t, wantQ012, uint192Big(q012))

	q0InvQ1, ok := inverseMod(q0%q1, q1)
	require.True(t, ok)
	q01Lo, q01Hi := bitsMul64ForTest(q0, q1)
	q01InvQ2, ok := inverseMod(mod128By64(q01Lo, q01Hi, q2), q2)
	require.True(t, ok)
	q012InvQ3, ok := inverseMod(mod192By64(q012, q3), q3)
	require.True(t, ok)

	values := []*big.Int{
		big.NewInt(0), big.NewInt(1), big.NewInt(-1),
		new(big.Int).Sub(new(big.Int).Set(wantHalf), big.NewInt(1)),
		new(big.Int).Set(wantHalf),
		new(big.Int).Add(new(big.Int).Set(wantHalf), big.NewInt(1)),
		new(big.Int).Neg(new(big.Int).Set(wantHalf)),
		new(big.Int).Add(new(big.Int).Neg(new(big.Int).Set(wantHalf)), big.NewInt(1)),
	}
	for i := int64(1); i <= 64; i++ {
		magnitude := new(big.Int).Lsh(big.NewInt(i*7919), uint((i*13)%150))
		values = append(values, magnitude, new(big.Int).Neg(new(big.Int).Set(magnitude)))
	}

	for _, signed := range values {
		var residues [MaxQPrefixWidth]uint64
		for row := 0; row < MaxQPrefixWidth; row++ {
			residues[row] = new(big.Int).Mod(new(big.Int).Set(signed), new(big.Int).SetUint64(q[row])).Uint64()
		}
		got := crtQ0123(residues[0], residues[1], residues[2], residues[3], q0, q1, q2, q3, q0InvQ1, q01InvQ2, q012InvQ3)
		wantCanonical := new(big.Int).Mod(new(big.Int).Set(signed), wantModulus)
		require.Equal(t, wantCanonical, uint192Big(got), "CRT value %s", signed)
		wantNegative := wantCanonical.Cmp(wantHalf) > 0
		wantMagnitude := new(big.Int).Set(wantCanonical)
		if wantNegative {
			wantMagnitude.Sub(wantModulus, wantMagnitude)
		}
		gotMagnitude, gotNegative := centeredQPrefix(got, modulus, half)
		require.Equal(t, wantNegative, gotNegative, "center sign %s", signed)
		require.Equal(t, wantMagnitude, uint192Big(gotMagnitude), "center magnitude %s", signed)
	}
}

func TestQ0123FixedWidthMultiplicationReportsOverflow(t *testing.T) {
	max := uint192{lo: ^uint64(0), mid: ^uint64(0), hi: ^uint64(0)}
	_, overflow := mul192By64(max, 3)
	require.True(t, overflow)
	_, _, _, overflow = q0123Modulus(^uint64(0), ^uint64(0)-2, ^uint64(0)-4, ^uint64(0)-6)
	require.True(t, overflow)
}
