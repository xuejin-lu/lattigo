package fast

import (
	"math/big"
	"testing"
)

// These benchmarks isolate the stages used by FastMulQ01. They intentionally
// use the same N=1024, level=8 fixture as the end-to-end multiplication
// benchmarks so their results can be compared directly.
func BenchmarkFastMulQ01Stages(b *testing.B) {
	r := benchmarkRing(b, 8)
	a := r.NewPoly()
	bb := r.NewPoly()
	for i := range a.Coeffs {
		for j := range a.Coeffs[i] {
			a.Coeffs[i][j] = uint64((j % 97) + 1)
			bb.Coeffs[i][j] = uint64((j % 53) + 1)
		}
	}

	aNTT := r.NewPoly()
	bNTT := r.NewPoly()
	cNTT := r.NewPoly()
	cCoeff := r.NewPoly()
	aMont := r.NewPoly()
	r.NTT(a, aNTT)
	r.NTT(bb, bNTT)
	r.SubRings[0].MForm(aNTT.Coeffs[0], aMont.Coeffs[0])
	r.SubRings[1].MForm(aNTT.Coeffs[1], aMont.Coeffs[1])
	r.SubRings[0].MulCoeffsMontgomery(aMont.Coeffs[0], bNTT.Coeffs[0], cNTT.Coeffs[0])
	r.SubRings[1].MulCoeffsMontgomery(aMont.Coeffs[1], bNTT.Coeffs[1], cNTT.Coeffs[1])
	if err := FastPartialINTT(r, cNTT, cCoeff); err != nil {
		b.Fatal(err)
	}

	bound := new(big.Int).Lsh(big.NewInt(1), 60)
	backend := NewRNSBackend(r, r.Level())

	b.Run("A_FastPartialNTT", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := FastPartialNTT(r, a, aNTT); err != nil {
				b.Fatal(err)
			}
		}
	})

	b.Run("B_MFormQ01", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			r.SubRings[0].MForm(aNTT.Coeffs[0], aMont.Coeffs[0])
			r.SubRings[1].MForm(aNTT.Coeffs[1], aMont.Coeffs[1])
		}
	})

	b.Run("C_MulCoeffsMontgomeryQ01", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			r.SubRings[0].MulCoeffsMontgomery(aMont.Coeffs[0], bNTT.Coeffs[0], cNTT.Coeffs[0])
			r.SubRings[1].MulCoeffsMontgomery(aMont.Coeffs[1], bNTT.Coeffs[1], cNTT.Coeffs[1])
		}
	})

	b.Run("D_FastPartialINTT", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := FastPartialINTT(r, cNTT, cCoeff); err != nil {
				b.Fatal(err)
			}
		}
	})

	b.Run("E_ReconstructQ0Q1", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, err := backend.ReconstructQ0Q1(&cCoeff); err != nil {
				b.Fatal(err)
			}
		}
	})

	coeffs, err := backend.ReconstructQ0Q1(&cCoeff)
	if err != nil {
		b.Fatal(err)
	}
	b.Run("F_CheckCoefficientBounds", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			for _, coeff := range coeffs {
				if err := backend.CheckCoefficientBounds(coeff, bound); err != nil {
					b.Fatal(err)
				}
			}
		}
	})

	b.Run("G_Redistribute", func(b *testing.B) {
		out := r.NewPoly()
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := backend.Redistribute(coeffs, &out); err != nil {
				b.Fatal(err)
			}
		}
	})
}
