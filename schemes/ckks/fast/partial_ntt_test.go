package fast

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func testPartialRing(t *testing.T, level int) *ring.Ring {
	t.Helper()
	const N = 16
	g := ring.NewNTTFriendlyPrimesGenerator(50, 2*N)
	moduli, err := g.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	r, err := ring.NewRing(N, moduli)
	require.NoError(t, err)
	return r
}

func fillPoly(p ring.Poly, seed int64) {
	rng := rand.New(rand.NewSource(seed))
	for i := range p.Coeffs {
		for j := range p.Coeffs[i] {
			p.Coeffs[i][j] = rng.Uint64() % 100000
		}
	}
}

func requireLimbEqual(t *testing.T, want, got []uint64) {
	t.Helper()
	require.Equal(t, want, got)
}

func TestFastPartialNTTMatchesSubRing(t *testing.T) {
	for _, level := range []int{1, 3} {
		r := testPartialRing(t, level)
		p := r.NewPoly()
		fillPoly(p, int64(level))
		want := r.NewPoly()
		r.SubRings[0].NTT(p.Coeffs[0], want.Coeffs[0])
		r.SubRings[1].NTT(p.Coeffs[1], want.Coeffs[1])
		for i := 2; i <= level; i++ {
			want.Coeffs[i] = make([]uint64, r.N())
			copy(want.Coeffs[i], p.Coeffs[i])
		}

		got := r.NewPoly()
		for i := 2; i <= level; i++ {
			got.Coeffs[i] = make([]uint64, r.N())
			copy(got.Coeffs[i], p.Coeffs[i])
		}
		require.NoError(t, FastPartialNTT(r, p, got))
		require.Equal(t, want.Coeffs[0], got.Coeffs[0])
		require.Equal(t, want.Coeffs[1], got.Coeffs[1])
		for i := 2; i <= level; i++ {
			requireLimbEqual(t, p.Coeffs[i], got.Coeffs[i])
		}
	}
}

func TestFastPartialINTTMatchesSubRing(t *testing.T) {
	for _, level := range []int{1, 3} {
		r := testPartialRing(t, level)
		p := r.NewPoly()
		fillPoly(p, int64(level+10))
		ntt := r.NewPoly()
		r.SubRings[0].NTT(p.Coeffs[0], ntt.Coeffs[0])
		r.SubRings[1].NTT(p.Coeffs[1], ntt.Coeffs[1])
		for i := 2; i <= level; i++ {
			copy(ntt.Coeffs[i], p.Coeffs[i])
		}

		want := r.NewPoly()
		r.SubRings[0].INTT(ntt.Coeffs[0], want.Coeffs[0])
		r.SubRings[1].INTT(ntt.Coeffs[1], want.Coeffs[1])
		got := r.NewPoly()
		for i := 2; i <= level; i++ {
			for j := range got.Coeffs[i] {
				got.Coeffs[i][j] = ^uint64(0) - uint64(i+j)
			}
		}
		before := make([][]uint64, level-1)
		for i := 2; i <= level; i++ {
			before[i-2] = append([]uint64(nil), got.Coeffs[i]...)
		}

		require.NoError(t, FastPartialINTT(r, ntt, got))
		require.Equal(t, want.Coeffs[0], got.Coeffs[0])
		require.Equal(t, want.Coeffs[1], got.Coeffs[1])
		for i := 2; i <= level; i++ {
			requireLimbEqual(t, before[i-2], got.Coeffs[i])
		}
	}
}

func TestFastPartialNTTRoundTripAndInPlace(t *testing.T) {
	r := testPartialRing(t, 3)
	p := r.NewPoly()
	fillPoly(p, 42)
	original := p.CopyNew()
	dormant := make([][]uint64, 2)
	for i := 2; i <= 3; i++ {
		dormant[i-2] = append([]uint64(nil), p.Coeffs[i]...)
	}

	require.NoError(t, FastPartialNTT(r, p, p))
	for i := 2; i <= 3; i++ {
		requireLimbEqual(t, dormant[i-2], p.Coeffs[i])
	}
	require.NoError(t, FastPartialINTT(r, p, p))
	require.Equal(t, original.Coeffs, p.Coeffs)
	for i := 2; i <= 3; i++ {
		requireLimbEqual(t, dormant[i-2], p.Coeffs[i])
	}
}

func TestFastPartialNTTEdgeCases(t *testing.T) {
	r := testPartialRing(t, 1)
	p := r.NewPoly()
	for i := range p.Coeffs[0] {
		var value int64
		switch i % 3 {
		case 0:
			value = 0
		case 1:
			value = 123
		default:
			value = -456
		}
		for limb := 0; limb < 2; limb++ {
			q := r.SubRings[limb].Modulus
			if value < 0 {
				p.Coeffs[limb][i] = q - uint64(-value)
			} else {
				p.Coeffs[limb][i] = uint64(value)
			}
		}
	}
	original := p.CopyNew()
	require.NoError(t, FastPartialNTT(r, p, p))
	require.NoError(t, FastPartialINTT(r, p, p))
	require.Equal(t, original.Coeffs, p.Coeffs)

	// INTT returns canonical coefficient residues suitable for the Phase 1A CRT.
	backend := NewRNSBackend(r, 1)
	coeffs, err := backend.ReconstructQ0Q1(&p)
	require.NoError(t, err)
	for i, coeff := range coeffs {
		switch i % 3 {
		case 0:
			require.Equal(t, big.NewInt(0), coeff)
		case 1:
			require.Equal(t, big.NewInt(123), coeff)
		default:
			require.Equal(t, big.NewInt(-456), coeff)
		}
	}
}

func BenchmarkFastPartialNTT(b *testing.B) {
	r := benchmarkRing(b, 8)
	p := r.NewPoly()
	fillPoly(p, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if err := FastPartialNTT(r, p, p); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkStandardFullLimbNTT(b *testing.B) {
	r := benchmarkRing(b, 8)
	p := r.NewPoly()
	fillPoly(p, 1)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.NTT(p, p)
	}
}

func BenchmarkFastPartialINTT(b *testing.B) {
	r := benchmarkRing(b, 8)
	p := r.NewPoly()
	fillPoly(p, 1)
	r.NTT(p, p)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if err := FastPartialINTT(r, p, p); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkStandardFullLimbINTT(b *testing.B) {
	r := benchmarkRing(b, 8)
	p := r.NewPoly()
	fillPoly(p, 1)
	r.NTT(p, p)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.INTT(p, p)
	}
}

func benchmarkRing(b *testing.B, level int) *ring.Ring {
	b.Helper()
	const N = 1024
	g := ring.NewNTTFriendlyPrimesGenerator(50, 2*N)
	moduli, err := g.NextAlternatingPrimes(level + 1)
	if err != nil {
		b.Fatal(err)
	}
	r, err := ring.NewRing(N, moduli)
	if err != nil {
		b.Fatal(err)
	}
	return r
}
