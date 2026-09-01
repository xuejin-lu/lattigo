package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func newFastTestCiphertext(t *testing.T, r *ring.Ring, degree, level int) *rlwe.Ciphertext {
	t.Helper()
	polys := make([]ring.Poly, degree+1)
	for i := range polys {
		polys[i] = ring.NewPoly(r.N(), level)
	}
	ct, err := rlwe.NewCiphertextAtLevelFromPoly(level, polys)
	require.NoError(t, err)
	ct.Scale = rlwe.NewScale(1)
	return ct
}

func fillFastCiphertext(ct *rlwe.Ciphertext, r *ring.Ring, offset uint64) {
	for d := range ct.Value {
		for i, s := range r.SubRings[:ct.Level()+1] {
			for j := range ct.Value[d].Coeffs[i] {
				ct.Value[d].Coeffs[i][j] = (offset + uint64(17*d+3*i+j)) % s.Modulus
			}
		}
	}
}

func TestFastAddSubQ01(t *testing.T) {
	for _, level := range []int{1, 2, 4} {
		r := testPartialRing(t, level)
		a := newFastTestCiphertext(t, r, 1, level)
		b := newFastTestCiphertext(t, r, 1, level)
		addOut := newFastTestCiphertext(t, r, 1, level)
		subOut := newFastTestCiphertext(t, r, 1, level)
		fillFastCiphertext(a, r, 11)
		fillFastCiphertext(b, r, 29)
		fillFastCiphertext(addOut, r, 101)
		fillFastCiphertext(subOut, r, 103)

		addDormant := cloneDormant(addOut, 2)
		subDormant := cloneDormant(subOut, 2)
		require.NoError(t, FastAdd(r, a, b, addOut))
		require.NoError(t, FastSub(r, a, b, subOut))

		for d := 0; d <= 1; d++ {
			for j := 0; j < r.N(); j++ {
				for _, i := range []int{0, 1} {
					q := r.SubRings[i].Modulus
					require.Equal(t, (a.Value[d].Coeffs[i][j]+b.Value[d].Coeffs[i][j])%q, addOut.Value[d].Coeffs[i][j])
					expectedSub := a.Value[d].Coeffs[i][j] - b.Value[d].Coeffs[i][j]
					if a.Value[d].Coeffs[i][j] < b.Value[d].Coeffs[i][j] {
						expectedSub += q
					}
					require.Equal(t, expectedSub, subOut.Value[d].Coeffs[i][j])
				}
			}
		}
		require.Equal(t, addDormant, cloneDormant(addOut, 2))
		require.Equal(t, subDormant, cloneDormant(subOut, 2))
	}
}

func TestFastAddSubIgnoresDormantInputs(t *testing.T) {
	r := testPartialRing(t, 3)
	a := newFastTestCiphertext(t, r, 1, 3)
	b := newFastTestCiphertext(t, r, 1, 3)
	a2 := a.CopyNew()
	b2 := b.CopyNew()
	fillFastCiphertext(a, r, 7)
	fillFastCiphertext(b, r, 19)
	*a2 = *a.CopyNew()
	*b2 = *b.CopyNew()
	for d := range a2.Value {
		for i := 2; i <= r.Level(); i++ {
			for j := range a2.Value[d].Coeffs[i] {
				a2.Value[d].Coeffs[i][j] = ^uint64(0) - uint64(i+j+d)
				b2.Value[d].Coeffs[i][j] = uint64(123 + i*j + j + d)
			}
		}
	}

	for _, sub := range []bool{false, true} {
		out1 := newFastTestCiphertext(t, r, 1, 3)
		out2 := newFastTestCiphertext(t, r, 1, 3)
		fillFastCiphertext(out1, r, 101)
		fillFastCiphertext(out2, r, 101)
		out1Dormant := cloneDormant(out1, 2)
		out2Dormant := cloneDormant(out2, 2)
		if sub {
			require.NoError(t, FastSub(r, a, b, out1))
			require.NoError(t, FastSub(r, a2, b2, out2))
		} else {
			require.NoError(t, FastAdd(r, a, b, out1))
			require.NoError(t, FastAdd(r, a2, b2, out2))
		}
		require.Equal(t, out1.Value[0].Coeffs[:2], out2.Value[0].Coeffs[:2])
		require.Equal(t, out1.Value[1].Coeffs[:2], out2.Value[1].Coeffs[:2])
		require.Equal(t, out1Dormant, cloneDormant(out1, 2))
		require.Equal(t, out2Dormant, cloneDormant(out2, 2))
	}
}

func TestFastAddSubInPlaceAndNTTDomain(t *testing.T) {
	r := testPartialRing(t, 2)
	a := newFastTestCiphertext(t, r, 1, 2)
	b := newFastTestCiphertext(t, r, 1, 2)
	fillFastCiphertext(a, r, 3)
	fillFastCiphertext(b, r, 41)

	want := b.CopyNew()
	require.NoError(t, FastAdd(r, a, b, want))
	require.NoError(t, FastAdd(r, a, b, a))
	require.Equal(t, want.Value[0].Coeffs[:2], a.Value[0].Coeffs[:2])
	require.Equal(t, want.Value[1].Coeffs[:2], a.Value[1].Coeffs[:2])

	left := a.CopyNew()
	right := b.CopyNew()
	wantSub := left.CopyNew()
	require.NoError(t, FastSub(r, left, right, wantSub))
	require.NoError(t, FastSub(r, left, right, right))
	require.Equal(t, wantSub.Value[0].Coeffs[:2], right.Value[0].Coeffs[:2])
	require.Equal(t, wantSub.Value[1].Coeffs[:2], right.Value[1].Coeffs[:2])

	nttA := newFastTestCiphertext(t, r, 1, 2)
	nttB := newFastTestCiphertext(t, r, 1, 2)
	fillFastCiphertext(nttA, r, 5)
	fillFastCiphertext(nttB, r, 13)
	for d := range nttA.Value {
		for i := 0; i < 2; i++ {
			r.SubRings[i].NTT(nttA.Value[d].Coeffs[i], nttA.Value[d].Coeffs[i])
			r.SubRings[i].NTT(nttB.Value[d].Coeffs[i], nttB.Value[d].Coeffs[i])
		}
	}
	nttA.IsNTT, nttB.IsNTT = true, true
	nttOut := newFastTestCiphertext(t, r, 1, 2)
	nttOut.IsNTT = true
	require.NoError(t, FastAdd(r, nttA, nttB, nttOut))
	require.True(t, nttOut.IsNTT)
}

func TestFastAddSubDifferentDegrees(t *testing.T) {
	r := testPartialRing(t, 2)
	low := newFastTestCiphertext(t, r, 1, 2)
	high := newFastTestCiphertext(t, r, 2, 2)
	out := newFastTestCiphertext(t, r, 2, 2)
	fillFastCiphertext(low, r, 3)
	fillFastCiphertext(high, r, 17)
	fillFastCiphertext(out, r, 101)
	outDormant := cloneDormant(out, 2)

	require.NoError(t, FastAdd(r, low, high, out))
	for j := 0; j < r.N(); j++ {
		for i := 0; i < 2; i++ {
			require.Equal(t, high.Value[2].Coeffs[i][j], out.Value[2].Coeffs[i][j])
		}
	}
	require.Equal(t, outDormant[0], cloneDormant(out, 2)[0])

	freshOut := newFastTestCiphertext(t, r, 2, 2)
	fillFastCiphertext(freshOut, r, 101)
	require.NoError(t, FastSub(r, low, high, freshOut))
	for j := 0; j < r.N(); j++ {
		for i := 0; i < 2; i++ {
			q := r.SubRings[i].Modulus
			expected := (q - high.Value[2].Coeffs[i][j]) % q
			require.Equal(t, expected, freshOut.Value[2].Coeffs[i][j])
		}
	}
}

func cloneDormant(ct *rlwe.Ciphertext, first int) [][][]uint64 {
	res := make([][][]uint64, len(ct.Value))
	for d := range ct.Value {
		res[d] = make([][]uint64, len(ct.Value[d].Coeffs)-first)
		for i := first; i < len(ct.Value[d].Coeffs); i++ {
			res[d][i-first] = append([]uint64(nil), ct.Value[d].Coeffs[i]...)
		}
	}
	return res
}

func BenchmarkFastAddSub(b *testing.B) {
	r := benchmarkRing(b, 8)
	a := newFastTestCiphertextBenchmark(r, 1, 8)
	bb := newFastTestCiphertextBenchmark(r, 1, 8)
	out := newFastTestCiphertextBenchmark(r, 1, 8)
	fillFastCiphertextBenchmark(a, r, 7)
	fillFastCiphertextBenchmark(bb, r, 19)
	b.Run("FastAdd", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := FastAdd(r, a, bb, out); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("FastSub", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := FastSub(r, a, bb, out); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("StandardAdd", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			for d := range a.Value {
				r.Add(a.Value[d], bb.Value[d], out.Value[d])
			}
		}
	})
	b.Run("StandardSub", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			for d := range a.Value {
				r.Sub(a.Value[d], bb.Value[d], out.Value[d])
			}
		}
	})
}

func newFastTestCiphertextBenchmark(r *ring.Ring, degree, level int) *rlwe.Ciphertext {
	polys := make([]ring.Poly, degree+1)
	for i := range polys {
		polys[i] = ring.NewPoly(r.N(), level)
	}
	ct, _ := rlwe.NewCiphertextAtLevelFromPoly(level, polys)
	ct.Scale = rlwe.NewScale(1)
	return ct
}

func fillFastCiphertextBenchmark(ct *rlwe.Ciphertext, r *ring.Ring, offset uint64) {
	for d := range ct.Value {
		for i, s := range r.SubRings[:ct.Level()+1] {
			for j := range ct.Value[d].Coeffs[i] {
				ct.Value[d].Coeffs[i][j] = (offset + uint64(j+3*i+17*d)) % s.Modulus
			}
		}
	}
}
