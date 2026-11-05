package ring

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
)

func compactTestRing(t *testing.T) *Ring {
	t.Helper()
	r, err := NewRing(16, Qi60[:6])
	require.NoError(t, err)
	return r
}

func requireCompactWidth(t *testing.T, p Poly, rows int) {
	t.Helper()
	for i, row := range p.Coeffs {
		if i < rows {
			require.Len(t, row, 16)
		} else {
			require.Empty(t, row)
		}
	}
}

func TestCompactRingModularOracleAndReuse(t *testing.T) {
	r := compactTestRing(t)
	for _, width := range []int{0, 1, 2, 4, 5, 6} {
		a, b := compactFixture(16, 5, width), compactFixture(16, 5, 6)
		out := compactFixture(16, 5, 6)
		r.Add(a, b, out)
		requireCompactWidth(t, out, width)
		for i := 0; i < width; i++ {
			for j, v := range out.Coeffs[i] {
				require.Equal(t, (a.Coeffs[i][j]+b.Coeffs[i][j])%r.SubRings[i].Modulus, v)
			}
		}
		// An overwrite must restore rows from complete inputs, unlike accumulation.
		r.Add(b, b, out)
		requireCompactWidth(t, out, 6)
		for i, row := range out.Coeffs {
			for j, v := range row {
				require.Equal(t, 2*b.Coeffs[i][j], v)
			}
		}
		before := a.CopyNew()
		shallow := a
		r.MulScalarBigint(a, big.NewInt(-7), shallow)
		requireCompactWidth(t, shallow, width)
		for i := 0; i < width; i++ {
			for j, v := range shallow.Coeffs[i] {
				want := new(big.Int).Mul(new(big.Int).SetUint64(before.Coeffs[i][j]), big.NewInt(-7))
				want.Mod(want, new(big.Int).SetUint64(r.SubRings[i].Modulus))
				require.Equal(t, want.Uint64(), v)
			}
		}
		// Independent modular oracle, rather than comparison with the same kernel.
		a = *before
		out.Zero()
		for repeat := 0; repeat < 3; repeat++ {
			r.MulCoeffsMontgomeryThenAdd(a, b, out)
			requireCompactWidth(t, out, width)
			for i := 0; i < width; i++ {
				q := new(big.Int).SetUint64(r.SubRings[i].Modulus)
				invR := new(big.Int).ModInverse(new(big.Int).Lsh(big.NewInt(1), 64), q)
				for j, v := range out.Coeffs[i] {
					want := new(big.Int).Mul(new(big.Int).SetUint64(a.Coeffs[i][j]), new(big.Int).SetUint64(b.Coeffs[i][j]))
					want.Mul(want, invR)
					want.Mul(want, big.NewInt(int64(repeat+1)))
					want.Mod(want, q)
					require.Equal(t, want.Uint64(), v)
				}
			}
		}
		out.Zero()
		r.Add(b, b, out)
		requireCompactWidth(t, out, 6)
	}
}

func TestCompactOverwriteAndAccumulatorAuthority(t *testing.T) {
	r := compactTestRing(t)
	a, b := compactFixture(16, 5, 4), compactFixture(16, 5, 6)
	mont := r.NewPoly()
	r.MForm(b, mont)
	out := compactFixture(16, 5, 1)
	r.MulCoeffsMontgomery(mont, a, out)
	requireCompactWidth(t, out, 4)
	for i, row := range out.Coeffs[:4] {
		for j, v := range row {
			require.Equal(t, b.Coeffs[i][j]*a.Coeffs[i][j], v)
		}
	}
	r.MulCoeffsMontgomery(mont, b, out)
	requireCompactWidth(t, out, 6)
	// Here the wider source itself is the accumulator/output via shallow alias.
	shallow := b
	savedHigh := b.Coeffs[5]
	before := b.CopyNew()
	r.MulCoeffsMontgomeryThenAdd(mont, a, shallow)
	requireCompactWidth(t, b, 4)
	for i, row := range b.Coeffs[:4] {
		for j, v := range row {
			require.Equal(t, before.Coeffs[i][j]+before.Coeffs[i][j]*a.Coeffs[i][j], v)
		}
	}
	// Dropping availability does not zero or mutate separately held backing.
	require.Equal(t, before.Coeffs[5], savedHigh)
}

func TestCompactOversizedAndPartialLevelRejectBeforeWriting(t *testing.T) {
	r := compactTestRing(t)
	a, out := compactFixture(16, 5, 4), compactFixture(32, 5, 6)
	before := out.CopyNew()
	require.Panics(t, func() { r.NTT(a, out) })
	require.Equal(t, *before, out)
	// An all-empty output must not be allocated too short for an oversized input.
	a = compactFixture(32, 5, 6)
	out = compactFixture(16, 5, 0)
	require.Panics(t, func() { r.MForm(a, out) })
	requireCompactWidth(t, out, 0)
	a = compactFixture(16, 5, 2)
	out = compactFixture(16, 5, 6)
	before = out.CopyNew()
	require.Panics(t, func() { out.CopyLvl(3, a) })
	require.Equal(t, *before, out)
	require.Panics(t, func() { r.AtLevel(3).Add(a, a, out) })
	require.Equal(t, *before, out)
}

func TestCompactRingDomainsPermutationAndAliases(t *testing.T) {
	r := compactTestRing(t)
	a := compactFixture(16, 5, 4)
	want := a.CopyNew()
	r.MForm(a, a)
	r.IMForm(a, a)
	require.True(t, a.Equal(want))
	for _, lazy := range []bool{false, true} {
		out := r.NewPoly()
		if lazy {
			r.NTTLazy(a, out)
			// NTTLazy may produce up to 6q; normalize before inverse input.
			for i := 0; i < 4; i++ {
				r.SubRings[i].Reduce(out.Coeffs[i], out.Coeffs[i])
			}
			r.INTTLazy(out, out)
		} else {
			r.NTT(a, out)
			r.INTT(out, out)
		}
		requireCompactWidth(t, out, 4)
		for i := 0; i < 4; i++ {
			for j, v := range out.Coeffs[i] {
				require.Equal(t, a.Coeffs[i][j], v%r.SubRings[i].Modulus)
			}
		}
	}
	index, err := AutomorphismNTTIndex(r.N(), r.NthRoot(), 5)
	require.NoError(t, err)
	out := r.NewPoly()
	r.AutomorphismNTTWithIndex(a, index, out)
	requireCompactWidth(t, out, 4)
	for i := 0; i < 4; i++ {
		for j, v := range out.Coeffs[i] {
			require.Equal(t, a.Coeffs[i][index[j]], v)
		}
	}
	// Exact coefficient backing alias remains valid for elementwise Add.
	view := Poly{Coeffs: append([][]uint64(nil), a.Coeffs...)}
	r.Add(a, a, view)
	for i := 0; i < 4; i++ {
		for j, v := range a.Coeffs[i] {
			require.Equal(t, 2*want.Coeffs[i][j], v)
		}
	}
	before := a.CopyNew()
	require.Panics(t, func() { r.AutomorphismNTTWithIndex(a, index, a) })
	require.True(t, a.Equal(before))
}

func TestCompactRingRejectionIsTransactional(t *testing.T) {
	r := compactTestRing(t)
	for _, kind := range []string{"hole", "short", "partial-alias", "cross-row", "headers"} {
		t.Run(kind, func(t *testing.T) {
			a, b, out := compactFixture(16, 5, 6), compactFixture(16, 5, 6), compactFixture(16, 5, 6)
			switch kind {
			case "hole":
				b.Coeffs[3] = nil
			case "short":
				b.Coeffs[4] = b.Coeffs[4][:15]
			case "partial-alias":
				backing := make([]uint64, 17)
				a.Coeffs[0] = backing[:16]
				out.Coeffs[0] = backing[1:]
			case "cross-row":
				out.Coeffs[3] = a.Coeffs[1]
			case "headers":
				headers := compactFixture(16, 6, 7).Coeffs
				a.Coeffs = headers[:6]
				out.Coeffs = headers[1:]
			}
			before := out.CopyNew()
			require.Panics(t, func() { r.Add(a, b, out) })
			require.Equal(t, *before, out)
		})
	}
	for _, kind := range []string{"length", "range", "duplicate", "alias"} {
		a, out := compactFixture(16, 5, 4), compactFixture(16, 5, 6)
		index, _ := AutomorphismNTTIndex(16, r.NthRoot(), 5)
		switch kind {
		case "length":
			index = index[:15]
		case "range":
			index[15] = 16
		case "duplicate":
			index[15] = index[0]
		case "alias":
			copy(out.Coeffs[5], index)
			index = out.Coeffs[5]
		}
		before := out.CopyNew()
		require.Panics(t, func() { r.AutomorphismNTTWithIndex(a, index, out) })
		require.Equal(t, *before, out)
	}
}
