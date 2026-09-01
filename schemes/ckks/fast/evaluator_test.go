package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func testFastCKKSParameters(t *testing.T) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{50, 50, 50, 50},
		LogDefaultScale: 40,
	})
	require.NoError(t, err)
	return params
}

func fillFastCKKSCiphertext(ct *rlwe.Ciphertext, params ckks.Parameters, values [][]int64, offset uint64) {
	r := params.RingQ().AtLevel(ct.Level())
	for d := range ct.Value {
		for i, s := range r.SubRings[:ct.Level()+1] {
			for j := range ct.Value[d].Coeffs[i] {
				if d < len(values) && j < len(values[d]) {
					v := values[d][j]
					if v < 0 {
						ct.Value[d].Coeffs[i][j] = s.Modulus - uint64(-v)
					} else {
						ct.Value[d].Coeffs[i][j] = uint64(v)
					}
				} else {
					ct.Value[d].Coeffs[i][j] = (offset + uint64(13*d+7*i+j)) % s.Modulus
				}
			}
		}
	}
}

func negacyclicProductMod(a, b []int64, q uint64) []uint64 {
	out := make([]uint64, len(a))
	for i := range a {
		var acc int64
		for j := range a {
			k := i - j
			sign := int64(1)
			if k < 0 {
				k += len(a)
				sign = -1
			}
			acc += sign * a[j] * b[k]
		}
		v := acc % int64(q)
		if v < 0 {
			v += int64(q)
		}
		out[i] = uint64(v)
	}
	return out
}

func TestFastMulQ01Authoritative(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	valuesA := [][]int64{{1, -2, 0, 3, -1, 2, 1, 0, -3, 2, 1, -1, 0, 2, -2, 1}, {2, 1, -1, 0, 3, -2, 1, 1, 0, -1, 2, 1, -2, 0, 1, 3}}
	valuesB := [][]int64{{2, 0, -1, 1, 1, -2, 3, 0, 1, 2, -1, 0, 2, -2, 1, 1}, {1, 2, 0, -1, 2, 1, -2, 1, 0, 3, -1, 2, 1, -2, 0, 1}}

	a := ckks.NewCiphertext(params, 1, 3)
	b := ckks.NewCiphertext(params, 1, 3)
	out := ckks.NewCiphertext(params, 2, 3)
	a.IsNTT, b.IsNTT, out.IsNTT = false, false, false
	fillFastCKKSCiphertext(a, params, valuesA, 17)
	fillFastCKKSCiphertext(b, params, valuesB, 31)
	fillFastCKKSCiphertext(out, params, nil, 101)
	dormant := cloneDormant(out, 2)

	require.NoError(t, eval.Mul(a, b, out))
	r := params.RingQ().AtLevel(3)
	for i := 0; i < 2; i++ {
		want0 := negacyclicProductMod(valuesA[0], valuesB[0], r.SubRings[i].Modulus)
		want2 := negacyclicProductMod(valuesA[1], valuesB[1], r.SubRings[i].Modulus)
		cross0 := negacyclicProductMod(valuesA[0], valuesB[1], r.SubRings[i].Modulus)
		cross1 := negacyclicProductMod(valuesA[1], valuesB[0], r.SubRings[i].Modulus)
		for j := 0; j < r.N(); j++ {
			require.Equal(t, want0[j], out.Value[0].Coeffs[i][j], "c0 limb=%d coeff=%d", i, j)
			require.Equal(t, (cross0[j]+cross1[j])%r.SubRings[i].Modulus, out.Value[1].Coeffs[i][j], "c1 limb=%d coeff=%d", i, j)
			require.Equal(t, want2[j], out.Value[2].Coeffs[i][j], "c2 limb=%d coeff=%d", i, j)
		}
	}
	require.Equal(t, 2, out.Degree())
	require.Equal(t, a.Scale.Mul(b.Scale), out.Scale)
	require.Equal(t, dormant, cloneDormant(out, 2))
}

func TestFastMulQ01AuthoritativeIgnoresDormantInputs(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	a := ckks.NewCiphertext(params, 1, 3)
	b := ckks.NewCiphertext(params, 1, 3)
	a2 := a.CopyNew()
	b2 := b.CopyNew()
	a.IsNTT, b.IsNTT, a2.IsNTT, b2.IsNTT = false, false, false, false
	fillFastCKKSCiphertext(a, params, nil, 7)
	fillFastCKKSCiphertext(b, params, nil, 19)
	*a2, *b2 = *a.CopyNew(), *b.CopyNew()
	for _, ct := range []*rlwe.Ciphertext{a2, b2} {
		for d := range ct.Value {
			for i := 2; i <= ct.Level(); i++ {
				for j := range ct.Value[d].Coeffs[i] {
					ct.Value[d].Coeffs[i][j] = ^uint64(0) - uint64(i+j+d)
				}
			}
		}
	}
	out1 := ckks.NewCiphertext(params, 2, 3)
	out2 := ckks.NewCiphertext(params, 2, 3)
	out1.IsNTT, out2.IsNTT = false, false
	fillFastCKKSCiphertext(out1, params, nil, 101)
	fillFastCKKSCiphertext(out2, params, nil, 103)
	dormant1 := cloneDormant(out1, 2)
	dormant2 := cloneDormant(out2, 2)
	require.NoError(t, eval.Mul(a, b, out1))
	require.NoError(t, eval.Mul(a2, b2, out2))
	for d := 0; d < 3; d++ {
		require.Equal(t, out1.Value[d].Coeffs[:2], out2.Value[d].Coeffs[:2])
	}
	require.Equal(t, dormant1, cloneDormant(out1, 2))
	require.Equal(t, dormant2, cloneDormant(out2, 2))
}

func TestFastEvaluatorMulInPlaceAndLevels(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	for _, level := range []int{1, 2, 3} {
		a := ckks.NewCiphertext(params, 1, level)
		b := ckks.NewCiphertext(params, 1, level)
		a.IsNTT, b.IsNTT = false, false
		fillFastCKKSCiphertext(a, params, nil, uint64(level+1))
		fillFastCKKSCiphertext(b, params, nil, uint64(level+11))
		want := a.CopyNew()
		require.NoError(t, eval.Mul(a, b, want))
		require.Equal(t, 2, want.Degree())
		require.Equal(t, level, want.Level())

		a2 := a.CopyNew()
		b2 := b.CopyNew()
		expected := ckks.NewCiphertext(params, 2, level)
		expected.IsNTT = false
		require.NoError(t, eval.Mul(a2, b2, expected))
		require.NoError(t, eval.Mul(a2, b2, a2))
		require.Equal(t, expected.Value[0].Coeffs[:2], a2.Value[0].Coeffs[:2])
		require.Equal(t, expected.Value[1].Coeffs[:2], a2.Value[1].Coeffs[:2])
		require.Equal(t, expected.Value[2].Coeffs[:2], a2.Value[2].Coeffs[:2])

		a3 := a.CopyNew()
		b3 := b.CopyNew()
		require.NoError(t, eval.Mul(a3, b3, b3))
		require.Equal(t, expected.Value[0].Coeffs[:2], b3.Value[0].Coeffs[:2])
		require.Equal(t, expected.Value[1].Coeffs[:2], b3.Value[1].Coeffs[:2])
		require.Equal(t, expected.Value[2].Coeffs[:2], b3.Value[2].Coeffs[:2])
	}
}

func BenchmarkFastEvaluatorMul(b *testing.B) {
	params := benchmarkFastCKKSParameters(b)
	fastEval := NewEvaluator(params)
	standardEval := ckks.NewEvaluator(params, nil)
	fastA := ckks.NewCiphertext(params, 1, params.MaxLevel())
	fastB := ckks.NewCiphertext(params, 1, params.MaxLevel())
	fastOut := ckks.NewCiphertext(params, 2, params.MaxLevel())
	fastA.IsNTT, fastB.IsNTT, fastOut.IsNTT = false, false, false
	fillFastCKKSCiphertextBenchmark(fastA, params, 7)
	fillFastCKKSCiphertextBenchmark(fastB, params, 19)
	standardA := fastA.CopyNew()
	standardB := fastB.CopyNew()
	standardOut := ckks.NewCiphertext(params, 2, params.MaxLevel())
	for d := range standardA.Value {
		params.RingQ().NTT(standardA.Value[d], standardA.Value[d])
		params.RingQ().NTT(standardB.Value[d], standardB.Value[d])
	}
	standardA.IsNTT, standardB.IsNTT, standardOut.IsNTT = true, true, true
	b.Run("FastAuthoritative", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := fastEval.Mul(fastA, fastB, fastOut); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("StandardEvaluator", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := standardEval.Mul(standardA, standardB, standardOut); err != nil {
				b.Fatal(err)
			}
		}
	})
}

func benchmarkFastCKKSParameters(b *testing.B) ckks.Parameters {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 10, LogQ: []int{50, 50, 50, 50, 50, 50, 50, 50, 50}, LogDefaultScale: 40})
	if err != nil {
		b.Fatal(err)
	}
	return params
}

func fillFastCKKSCiphertextBenchmark(ct *rlwe.Ciphertext, params ckks.Parameters, offset uint64) {
	for d := range ct.Value {
		for i, s := range params.RingQ().SubRings[:ct.Level()+1] {
			for j := range ct.Value[d].Coeffs[i] {
				ct.Value[d].Coeffs[i][j] = (offset + uint64(j+3*i+17*d)) % s.Modulus
			}
		}
	}
}
