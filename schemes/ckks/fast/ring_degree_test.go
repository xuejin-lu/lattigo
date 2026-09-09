package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func testFastRingPair(t *testing.T, level int) (small, large *ring.Ring) {
	t.Helper()
	const nSmall = 16
	generator := ring.NewNTTFriendlyPrimesGenerator(50, 2*2*nSmall)
	moduli, err := generator.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	small, err = ring.NewRing(nSmall, moduli)
	require.NoError(t, err)
	large, err = ring.NewRing(2*nSmall, moduli)
	require.NoError(t, err)
	return
}

func TestFastN1ToN2CoefficientDomain(t *testing.T) {
	for _, level := range []int{1, 2, 4} {
		n1, n2 := testFastRingPair(t, level)
		in := newFastTestCiphertext(t, n1, 1, level)
		out := newFastTestCiphertext(t, n2, 1, level)
		fillFastCiphertext(in, n1, 7)
		fillFastCiphertext(out, n2, 101)
		in.IsNTT, out.IsNTT = false, false
		in.IsMontgomery, out.IsMontgomery = true, true
		in.IsBatched, in.IsBitReversed = true, true
		in.LogDimensions = ring.Dimensions{Rows: 1, Cols: 3}
		in.Scale = rlwe.NewScale(123.5)
		metadata := *in.MetaData
		dormant := cloneDormant(out, 2)

		want := newFastTestCiphertext(t, n2, 1, level)
		rlwe.SwitchCiphertextRingDegree(in.El(), want.El())
		require.NoError(t, FastN1ToN2(n1, n2, in, out))
		require.Equal(t, want.Value[0].Coeffs[:2], out.Value[0].Coeffs[:2])
		require.Equal(t, want.Value[1].Coeffs[:2], out.Value[1].Coeffs[:2])
		require.Equal(t, metadata, *out.MetaData)
		require.Equal(t, dormant, cloneDormant(out, 2))
	}
}

func TestFastRingDegreeLevelZero(t *testing.T) {
	n1, n2 := testFastRingPair(t, 0)
	for _, expand := range []bool{true, false} {
		inRing, outRing := n1, n2
		if !expand {
			inRing, outRing = n2, n1
		}
		in := newFastTestCiphertext(t, inRing, 1, 0)
		out := newFastTestCiphertext(t, outRing, 1, 0)
		fillFastCiphertext(in, inRing, 17)
		in.IsNTT, out.IsNTT = false, false
		in.IsMontgomery, out.IsMontgomery = false, false
		if expand {
			require.NoError(t, FastN1ToN2(n1, n2, in, out))
		} else {
			require.NoError(t, FastN2ToN1(n2, n1, in, out))
		}
		for d := 0; d <= 1; d++ {
			if expand {
				for i, value := range in.Value[d].Coeffs[0] {
					require.Equal(t, value, out.Value[d].Coeffs[0][2*i])
					require.Zero(t, out.Value[d].Coeffs[0][2*i+1])
				}
			} else {
				for i, value := range out.Value[d].Coeffs[0] {
					require.Equal(t, value, in.Value[d].Coeffs[0][2*i])
				}
			}
		}
	}

	in := newFastTestCiphertext(t, n1, 1, 0)
	out := newFastTestCiphertext(t, n2, 1, 0)
	fillFastCiphertext(in, n1, 23)
	in.IsNTT, out.IsNTT = true, true
	in.IsMontgomery, out.IsMontgomery = true, true
	n1.SubRings[0].NTT(in.Value[0].Coeffs[0], in.Value[0].Coeffs[0])
	n1.SubRings[0].NTT(in.Value[1].Coeffs[0], in.Value[1].Coeffs[0])
	require.NoError(t, FastN1ToN2(n1, n2, in, out))
	require.True(t, out.IsNTT)
	require.True(t, out.IsMontgomery)
}

func TestFastN2ToN1CoefficientDomain(t *testing.T) {
	for _, level := range []int{1, 2, 4} {
		n1, n2 := testFastRingPair(t, level)
		in := newFastTestCiphertext(t, n2, 1, level)
		out := newFastTestCiphertext(t, n1, 1, level)
		fillFastCiphertext(in, n2, 17)
		fillFastCiphertext(out, n1, 201)
		in.IsNTT, out.IsNTT = false, false
		in.IsMontgomery, out.IsMontgomery = false, false
		in.IsBatched, in.IsBitReversed = true, true
		in.LogDimensions = ring.Dimensions{Rows: 1, Cols: 3}
		in.Scale = rlwe.NewScale(123.5)
		metadata := *in.MetaData
		dormant := cloneDormant(out, 2)

		want := newFastTestCiphertext(t, n1, 1, level)
		rlwe.SwitchCiphertextRingDegree(in.El(), want.El())
		require.NoError(t, FastN2ToN1(n2, n1, in, out))
		require.Equal(t, want.Value[0].Coeffs[:2], out.Value[0].Coeffs[:2])
		require.Equal(t, want.Value[1].Coeffs[:2], out.Value[1].Coeffs[:2])
		require.Equal(t, metadata, *out.MetaData)
		require.Equal(t, dormant, cloneDormant(out, 2))
	}
}

func TestFastRingDegreeNTTDomain(t *testing.T) {
	for _, expand := range []bool{true, false} {
		for _, level := range []int{1, 2, 4} {
			n1, n2 := testFastRingPair(t, level)
			inRing, outRing := n1, n2
			if !expand {
				inRing, outRing = n2, n1
			}
			in := newFastTestCiphertext(t, inRing, 1, level)
			out := newFastTestCiphertext(t, outRing, 1, level)
			fillFastCiphertext(in, inRing, 31)
			fillFastCiphertext(out, outRing, 401)
			in.IsNTT, out.IsNTT = true, true
			in.IsMontgomery, out.IsMontgomery = true, true
			for d := 0; d <= 1; d++ {
				for limb := 0; limb < 2; limb++ {
					inRing.SubRings[limb].NTT(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
				}
			}
			dormant := cloneDormant(out, 2)

			want := newFastTestCiphertext(t, outRing, 1, level)
			var largeRing *ring.Ring
			if expand {
				largeRing = nil
			} else {
				largeRing = n2
			}
			rlwe.SwitchCiphertextRingDegreeNTT(in.El(), largeRing, want.El())
			var err error
			if expand {
				err = FastN1ToN2(n1, n2, in, out)
			} else {
				err = FastN2ToN1(n2, n1, in, out)
			}
			require.NoError(t, err)
			require.Equal(t, want.Value[0].Coeffs[:2], out.Value[0].Coeffs[:2])
			require.Equal(t, want.Value[1].Coeffs[:2], out.Value[1].Coeffs[:2])
			require.True(t, out.IsNTT)
			require.True(t, out.IsMontgomery)
			require.Equal(t, dormant, cloneDormant(out, 2))
		}
	}
}

func TestFastRingDegreeIgnoresDormantLimbs(t *testing.T) {
	for _, expand := range []bool{true, false} {
		n1, n2 := testFastRingPair(t, 4)
		inRing, outRing := n1, n2
		if !expand {
			inRing, outRing = n2, n1
		}
		inA := newFastTestCiphertext(t, inRing, 1, 4)
		inB := inA.CopyNew()
		fillFastCiphertext(inA, inRing, 11)
		*inB = *inA.CopyNew()
		for d := range inB.Value {
			for limb := 2; limb <= inB.Level(); limb++ {
				for j := range inB.Value[d].Coeffs[limb] {
					inB.Value[d].Coeffs[limb][j] = uint64(10000 + d + limb + j)
				}
			}
		}
		outA := newFastTestCiphertext(t, outRing, 1, 4)
		outB := newFastTestCiphertext(t, outRing, 1, 4)
		fillFastCiphertext(outA, outRing, 71)
		fillFastCiphertext(outB, outRing, 71)
		dormantA, dormantB := cloneDormant(outA, 2), cloneDormant(outB, 2)
		var errA, errB error
		if expand {
			errA = FastN1ToN2(n1, n2, inA, outA)
			errB = FastN1ToN2(n1, n2, inB, outB)
		} else {
			errA = FastN2ToN1(n2, n1, inA, outA)
			errB = FastN2ToN1(n2, n1, inB, outB)
		}
		require.NoError(t, errA)
		require.NoError(t, errB)
		require.Equal(t, outA.Value[0].Coeffs[:2], outB.Value[0].Coeffs[:2])
		require.Equal(t, outA.Value[1].Coeffs[:2], outB.Value[1].Coeffs[:2])
		require.Equal(t, dormantA, cloneDormant(outA, 2))
		require.Equal(t, dormantB, cloneDormant(outB, 2))
	}
}

func TestFastRingDegreeRejectsAliasAndInvalidRings(t *testing.T) {
	n1, n2 := testFastRingPair(t, 1)
	ct1 := newFastTestCiphertext(t, n1, 1, 1)
	ct2 := newFastTestCiphertext(t, n2, 1, 1)
	require.Error(t, FastN1ToN2(n1, n2, ct1, ct1))
	require.Error(t, FastN2ToN1(n2, n1, ct2, ct2))
	require.Error(t, FastN1ToN2(n2, n1, ct1, ct2))
	require.Error(t, FastN1ToN2(n1.AtLevel(0), n2.AtLevel(0), ct1, ct2))
}

func BenchmarkFastRingDegreeConversion(b *testing.B) {
	n1, n2 := benchmarkFastRingPair(b)
	in1 := newFastTestCiphertextBenchmark(n1, 1, 8)
	out2 := newFastTestCiphertextBenchmark(n2, 1, 8)
	in2 := newFastTestCiphertextBenchmark(n2, 1, 8)
	out1 := newFastTestCiphertextBenchmark(n1, 1, 8)
	standardOut2 := newFastTestCiphertextBenchmark(n2, 1, 8)
	standardOut1 := newFastTestCiphertextBenchmark(n1, 1, 8)
	fillFastCiphertextBenchmark(in1, n1, 7)
	fillFastCiphertextBenchmark(in2, n2, 19)

	b.Run("FastN1ToN2Coefficient", func(b *testing.B) {
		in1.IsNTT, out2.IsNTT = false, false
		for i := 0; i < b.N; i++ {
			if err := FastN1ToN2(n1, n2, in1, out2); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("StandardN1ToN2Coefficient", func(b *testing.B) {
		in1.IsNTT, standardOut2.IsNTT = false, false
		for i := 0; i < b.N; i++ {
			rlwe.SwitchCiphertextRingDegree(in1.El(), standardOut2.El())
		}
	})
	b.Run("FastN2ToN1Coefficient", func(b *testing.B) {
		in2.IsNTT, out1.IsNTT = false, false
		for i := 0; i < b.N; i++ {
			if err := FastN2ToN1(n2, n1, in2, out1); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("StandardN2ToN1Coefficient", func(b *testing.B) {
		in2.IsNTT, standardOut1.IsNTT = false, false
		for i := 0; i < b.N; i++ {
			rlwe.SwitchCiphertextRingDegree(in2.El(), standardOut1.El())
		}
	})

	nttIn1 := in1.CopyNew()
	nttIn2 := in2.CopyNew()
	nttOut2 := newFastTestCiphertextBenchmark(n2, 1, 8)
	nttOut1 := newFastTestCiphertextBenchmark(n1, 1, 8)
	for _, item := range []struct {
		ct *rlwe.Ciphertext
		r  *ring.Ring
	}{
		{nttIn1, n1},
		{nttIn2, n2},
	} {
		for d := 0; d <= 1; d++ {
			for limb := 0; limb < 2; limb++ {
				item.r.SubRings[limb].NTT(item.ct.Value[d].Coeffs[limb], item.ct.Value[d].Coeffs[limb])
			}
		}
		item.ct.IsNTT = true
	}
	nttOut2.IsNTT, nttOut1.IsNTT = true, true
	b.Run("FastN1ToN2NTT", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err := FastN1ToN2(n1, n2, nttIn1, nttOut2); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("StandardN1ToN2NTT", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rlwe.SwitchCiphertextRingDegreeNTT(nttIn1.El(), nil, nttOut2.El())
		}
	})
	b.Run("FastN2ToN1NTT", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err := FastN2ToN1(n2, n1, nttIn2, nttOut1); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("StandardN2ToN1NTT", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			rlwe.SwitchCiphertextRingDegreeNTT(nttIn2.El(), n2, nttOut1.El())
		}
	})
}

func benchmarkFastRingPair(b *testing.B) (small, large *ring.Ring) {
	b.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(50, 2*32)
	moduli, err := generator.NextAlternatingPrimes(9)
	if err != nil {
		b.Fatal(err)
	}
	small, err = ring.NewRing(16, moduli)
	if err != nil {
		b.Fatal(err)
	}
	large, err = ring.NewRing(32, moduli)
	if err != nil {
		b.Fatal(err)
	}
	return
}
