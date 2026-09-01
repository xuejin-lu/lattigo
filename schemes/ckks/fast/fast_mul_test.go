package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func smallCoefficientPoly(t *testing.T, r *ring.Ring, values []int64) ring.Poly {
	t.Helper()
	p := r.NewPoly()
	for k, value := range values {
		for i, subring := range r.SubRings[:r.Level()+1] {
			if value < 0 {
				p.Coeffs[i][k] = subring.Modulus - uint64(-value)
			} else {
				p.Coeffs[i][k] = uint64(value)
			}
		}
	}
	return p
}

func standardMul(t *testing.T, r *ring.Ring, p1, p2 ring.Poly) ring.Poly {
	t.Helper()
	aNTT := r.NewPoly()
	bNTT := r.NewPoly()
	cNTT := r.NewPoly()
	c := r.NewPoly()
	r.NTT(p1, aNTT)
	r.NTT(p2, bNTT)
	r.MForm(aNTT, aNTT)
	r.MulCoeffsMontgomery(aNTT, bNTT, cNTT)
	r.INTT(cNTT, c)
	return c
}

func requirePolyEqual(t *testing.T, want, got ring.Poly) {
	t.Helper()
	require.Equal(t, want.Coeffs, got.Coeffs)
}

func TestFastMulQ01MatchesStandard(t *testing.T) {
	for _, level := range []int{1, 2, 4} {
		r := testPartialRing(t, level)
		valuesA := make([]int64, r.N())
		valuesB := make([]int64, r.N())
		for i := range valuesA {
			valuesA[i] = int64((i % 7) - 3)
			valuesB[i] = int64((i % 5) - 2)
		}
		a := smallCoefficientPoly(t, r, valuesA)
		b := smallCoefficientPoly(t, r, valuesB)
		got := r.NewPoly()
		require.NoError(t, FastMulQ01(r, a, b, got, big.NewInt(1<<20)))
		want := standardMul(t, r, a, b)
		requirePolyEqual(t, want, got)

		backend := NewRNSBackend(r, level)
		gotCoeffs, err := backend.ReconstructQ0Q1(&got)
		require.NoError(t, err)
		wantCoeffs, err := backend.ReconstructQ0Q1(&want)
		require.NoError(t, err)
		for i := range gotCoeffs {
			require.Equal(t, wantCoeffs[i].String(), gotCoeffs[i].String())
		}
	}
}

func TestFastMulQ01IgnoresDormantInputLimbs(t *testing.T) {
	r := testPartialRing(t, 3)
	a := smallCoefficientPoly(t, r, []int64{1, -2, 0, 3, -1, 2, 1, 0, -3, 2, 1, -1, 0, 2, -2, 1})
	b := smallCoefficientPoly(t, r, []int64{2, 1, -1, 0, 3, -2, 1, 1, 0, -1, 2, 1, -2, 0, 1, 3})
	a2 := a.CopyNew()
	b2 := b.CopyNew()
	for i := 2; i <= r.Level(); i++ {
		for j := range a2.Coeffs[i] {
			a2.Coeffs[i][j] = ^uint64(0) - uint64(i+j)
			b2.Coeffs[i][j] = uint64(i*j + j + 17)
		}
	}
	aDormant := make([][]uint64, r.Level()-1)
	bDormant := make([][]uint64, r.Level()-1)
	for i := 2; i <= r.Level(); i++ {
		aDormant[i-2] = append([]uint64(nil), a.Coeffs[i]...)
		bDormant[i-2] = append([]uint64(nil), b.Coeffs[i]...)
	}

	got1 := r.NewPoly()
	got2 := r.NewPoly()
	require.NoError(t, FastMulQ01(r, a, b, got1, big.NewInt(1<<20)))
	require.NoError(t, FastMulQ01(r, *a2, *b2, got2, big.NewInt(1<<20)))
	requirePolyEqual(t, got1, got2)
	for i := 2; i <= r.Level(); i++ {
		require.Equal(t, aDormant[i-2], a.Coeffs[i])
		require.Equal(t, bDormant[i-2], b.Coeffs[i])
	}
}

func TestFastMulQ01OutputRedistribution(t *testing.T) {
	r := testPartialRing(t, 2)
	a := smallCoefficientPoly(t, r, []int64{1, 2, -1, 0, 3, -2, 1, 2, 0, -1, 2, 1, -2, 0, 1, 3})
	b := smallCoefficientPoly(t, r, []int64{2, 0, -1, 1, 1, -2, 3, 0, 1, 2, -1, 0, 2, -2, 1, 1})
	got := r.NewPoly()
	require.NoError(t, FastMulQ01(r, a, b, got, big.NewInt(1<<20)))
	backend := NewRNSBackend(r, r.Level())
	coeffs, err := backend.ReconstructQ0Q1(&got)
	require.NoError(t, err)
	for i, coeff := range coeffs {
		for limb := 0; limb <= r.Level(); limb++ {
			modulus := new(big.Int).SetUint64(r.SubRings[limb].Modulus)
			want := new(big.Int).Mod(coeff, modulus)
			require.Equal(t, want.Uint64(), got.Coeffs[limb][i], "limb=%d coeff=%d", limb, i)
		}
	}
}

func TestFastMulQ01RejectsInvalidInputs(t *testing.T) {
	r := testPartialRing(t, 1)
	p := r.NewPoly()
	require.Error(t, FastMulQ01(r, p, p, p, nil))
	require.Error(t, FastMulQ01(r, p, p, p, big.NewInt(-1)))

	r0, err := ring.NewRing(r.N(), []uint64{r.SubRings[0].Modulus})
	require.NoError(t, err)
	require.Error(t, FastMulQ01(r0, r0.NewPoly(), r0.NewPoly(), r0.NewPoly(), big.NewInt(1)))
}

func BenchmarkFastMulQ01(b *testing.B) {
	r := benchmarkRing(b, 8)
	a := r.NewPoly()
	bb := r.NewPoly()
	for i := range a.Coeffs {
		for j := range a.Coeffs[i] {
			a.Coeffs[i][j] = uint64((j % 97) + 1)
			bb.Coeffs[i][j] = uint64((j % 53) + 1)
		}
	}
	out := r.NewPoly()
	bound := new(big.Int).Lsh(big.NewInt(1), 45)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if err := FastMulQ01(r, a, bb, out, bound); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkStandardPolynomialMul(b *testing.B) {
	r := benchmarkRing(b, 8)
	a := r.NewPoly()
	bb := r.NewPoly()
	aNTT := r.NewPoly()
	bNTT := r.NewPoly()
	aMont := r.NewPoly()
	cNTT := r.NewPoly()
	out := r.NewPoly()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.NTT(a, aNTT)
		r.NTT(bb, bNTT)
		r.MForm(aNTT, aMont)
		r.MulCoeffsMontgomery(aMont, bNTT, cNTT)
		r.INTT(cNTT, out)
	}
}
