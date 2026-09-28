package fast

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestQPrefixFixedWidthDivisionsMatchBigInt(t *testing.T) {
	params := q012TestParameters(t)
	moduli := append([]uint64(nil), params.Q()...)
	moduli = append(moduli,
		1, 2, 3,
		(uint64(1)<<39)-1,
		(uint64(1)<<40)-87,
		(uint64(1)<<55)-19,
		(uint64(1)<<56)-5,
		(uint64(1)<<60)-1,
		(uint64(1)<<61)-1,
		^uint64(0),
	)
	rng := rand.New(rand.NewSource(0x51505245464958))
	max192 := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), 192), big.NewInt(1))

	for _, modulus := range moduli {
		t.Run(big.NewInt(0).SetUint64(modulus).String(), func(t *testing.T) {
			m := new(big.Int).SetUint64(modulus)
			values := []uint192{
				{}, {lo: 1}, {mid: 1}, {hi: 1},
				{lo: ^uint64(0)}, {mid: ^uint64(0)}, {hi: ^uint64(0)},
				{lo: ^uint64(0), mid: ^uint64(0)},
				{lo: ^uint64(0), mid: ^uint64(0), hi: ^uint64(0)},
			}
			for _, multiple := range []uint64{1, 2, 17, 4093} {
				base := new(big.Int).Mul(new(big.Int).Set(m), new(big.Int).SetUint64(multiple))
				for _, delta := range []int64{-1, 0, 1} {
					candidate := new(big.Int).Add(new(big.Int).Set(base), big.NewInt(delta))
					if candidate.Sign() >= 0 && candidate.Cmp(max192) <= 0 {
						values = append(values, bigUint192(candidate))
					}
				}
			}
			for i := 0; i < 128; i++ {
				values = append(values, uint192{lo: rng.Uint64(), mid: rng.Uint64(), hi: rng.Uint64()})
				values = append(values, uint192{lo: rng.Uint64(), mid: rng.Uint64()})
			}
			for _, value := range values {
				want := new(big.Int).Mod(uint192Big(value), m).Uint64()
				require.Equal(t, want, mod192By64(value, modulus), "value=%v modulus=%d", value, modulus)
				if value.hi == 0 {
					want128 := new(big.Int).Mod(uint192Big(value), m).Uint64()
					require.Equal(t, want128, mod128By64(value.lo, value.mid, modulus), "value=%v modulus=%d", value, modulus)
				}
			}
		})
	}
}

func TestRoundedMagnitude192MatchesBigIntAndExactDiv64Candidate(t *testing.T) {
	params := q012TestParameters(t)
	divisors := append([]uint64(nil), params.Q()...)
	divisors = append(divisors, 1, 2, 3, uint64(1)<<31, uint64(1)<<32, uint64(1)<<39, uint64(1)<<40, uint64(1)<<55, uint64(1)<<56, uint64(1)<<60, ^uint64(0))
	rng := rand.New(rand.NewSource(0x524f554e444d4147))
	max192 := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), 192), big.NewInt(1))

	for _, divisor := range divisors {
		t.Run(big.NewInt(0).SetUint64(divisor).String(), func(t *testing.T) {
			d := new(big.Int).SetUint64(divisor)
			values := []*big.Int{big.NewInt(0), big.NewInt(1), new(big.Int).Sub(new(big.Int).Set(d), big.NewInt(1)), new(big.Int).Set(d)}
			for _, multiplier := range []*big.Int{big.NewInt(0), big.NewInt(1), new(big.Int).Lsh(big.NewInt(1), 64), new(big.Int).Lsh(big.NewInt(1), 128)} {
				base := new(big.Int).Mul(new(big.Int).Set(multiplier), d)
				for _, delta := range []int64{-1, 0, 1} {
					center := new(big.Int).Add(new(big.Int).Set(base), new(big.Int).Rsh(new(big.Int).Set(d), 1))
					candidate := new(big.Int).Add(center, big.NewInt(delta))
					if candidate.Sign() >= 0 && candidate.Cmp(max192) <= 0 {
						values = append(values, candidate)
					}
				}
			}
			for i := 0; i < 128; i++ {
				values = append(values, uint192Big(uint192{lo: rng.Uint64(), mid: rng.Uint64(), hi: rng.Uint64()}))
			}
			for _, value := range values {
				input := bigUint192(value)
				want, remainder := new(big.Int).QuoRem(new(big.Int).Set(value), d, new(big.Int))
				if new(big.Int).Lsh(remainder, 1).Cmp(d) > 0 {
					want.Add(want, big.NewInt(1))
				}
				got := roundedMagnitude192(input, divisor)
				require.Zero(t, want.Cmp(uint192Big(got)), "value=%s divisor=%d", value, divisor)
				require.Equal(t, got, roundedMagnitude192BitsDiv64Candidate(input, divisor), "Div64 candidate value=%s divisor=%d", value, divisor)
			}
		})
	}
}

func TestReconstructQPrefixMatchesGenericCRTAndBigInt(t *testing.T) {
	params := q012TestParameters(t)
	q := params.Q()
	scratch := newFastRescaleScratch(params.RingQ())
	rng := rand.New(rand.NewSource(0x4352545052454649))
	for rows := 1; rows <= MaxQPrefixWidth; rows++ {
		product := big.NewInt(1)
		for _, modulus := range q[:rows] {
			product.Mul(product, new(big.Int).SetUint64(modulus))
		}
		values := []*big.Int{
			big.NewInt(0), big.NewInt(1), new(big.Int).Rsh(new(big.Int).Set(product), 1),
			new(big.Int).Sub(new(big.Int).Set(product), big.NewInt(1)),
		}
		for i := 0; i < 128; i++ {
			randomValue := uint192{lo: rng.Uint64(), mid: rng.Uint64(), hi: rng.Uint64()}
			values = append(values, new(big.Int).Mod(uint192Big(randomValue), product))
		}
		for _, value := range values {
			var residues [MaxQPrefixWidth]uint64
			for row := 0; row < rows; row++ {
				residues[row] = new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(q[row])).Uint64()
			}
			got := reconstructQPrefix(rows, residues, &scratch)
			var generic uint192
			switch rows {
			case 1:
				generic = uint192{lo: residues[0]}
			case 2:
				lo, hi := crtQ01(residues[0], residues[1], q[0], q[1], scratch.inverse[1])
				generic = uint192{lo: lo, mid: hi}
			case 3:
				generic = crtQ012(residues[0], residues[1], residues[2], q[0], q[1], q[2], scratch.inverse[1], scratch.inverse[2])
			case 4:
				generic = crtQ0123(residues[0], residues[1], residues[2], residues[3], q[0], q[1], q[2], q[3], scratch.inverse[1], scratch.inverse[2], scratch.inverse[3])
			}
			require.Equal(t, generic, got, "rows=%d value=%s", rows, value)
			require.Zero(t, value.Cmp(uint192Big(got)), "BigInt rows=%d", rows)
		}
	}
}
