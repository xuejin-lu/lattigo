package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func traceTestParameters(t *testing.T) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{50, 50, 50, 50, 50},
		LogDefaultScale: 40,
	})
	require.NoError(t, err)
	return params
}

func newTraceCiphertext(params ckks.Parameters, level int, montgomery, poison bool) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, level)
	ct.IsNTT = true
	ct.IsMontgomery = montgomery
	ct.Scale = rlwe.NewScale(123.5)
	for component := range ct.Value {
		for limb := 0; limb < 2 && limb <= level; limb++ {
			q := params.RingQ().SubRings[limb].Modulus
			for j := range ct.Value[component].Coeffs[limb] {
				ct.Value[component].Coeffs[limb][j] = (uint64(17+component*31+limb*43) + uint64(j*j+7*j)) % q
			}
			params.RingQ().SubRings[limb].NTT(ct.Value[component].Coeffs[limb], ct.Value[component].Coeffs[limb])
			if montgomery {
				params.RingQ().SubRings[limb].MForm(ct.Value[component].Coeffs[limb], ct.Value[component].Coeffs[limb])
			}
		}
		for limb := 2; limb <= level; limb++ {
			for j := range ct.Value[component].Coeffs[limb] {
				if poison {
					ct.Value[component].Coeffs[limb][j] = uint64(0x600d0000 + component + limb + j)
				} else {
					ct.Value[component].Coeffs[limb][j] = uint64(1000 + component + limb + j)
				}
			}
		}
	}
	return ct
}

func traceReference(params ckks.Parameters, ctIn *rlwe.Ciphertext, logN int) *rlwe.Ciphertext {
	got := ctIn.CopyNew()
	level := ctIn.Level()
	ringQ := params.RingQ().AtLevel(level)
	gap := 1 << (params.LogN() - logN - 1)
	if logN == 0 {
		gap <<= 1
	}
	if gap <= 1 {
		return got
	}
	nInv := new(big.Int).SetUint64(uint64(gap))
	requireTraceInverse(nInv, ringQ.ModulusAtLevel[level])
	for component := range got.Value {
		ringQ.MulScalarBigint(got.Value[component], nInv, got.Value[component])
	}
	apply := func(galEl uint64) {
		index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galEl)
		if err != nil {
			panic(err)
		}
		r := params.RingQ().AtLevel(1)
		for component := range got.Value {
			tmp := ring.NewPoly(r.N(), 1)
			r.AutomorphismNTTWithIndex(got.Value[component], index, tmp)
			for limb := 0; limb < 2; limb++ {
				r.SubRings[limb].Add(got.Value[component].Coeffs[limb], tmp.Coeffs[limb], got.Value[component].Coeffs[limb])
			}
		}
	}
	for i := logN; i < params.LogN()-1; i++ {
		apply(params.GaloisElement(1 << i))
	}
	if logN == 0 {
		apply(ringQ.NthRoot() - 1)
	}
	return got
}

func requireTraceInverse(x, modulus *big.Int) {
	if x.ModInverse(x, modulus) == nil {
		panic("trace test inverse does not exist")
	}
}

func requireTraceMatches(t *testing.T, want, got *rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, want.Level(), got.Level())
	require.Equal(t, want.Scale, got.Scale)
	require.Equal(t, want.IsNTT, got.IsNTT)
	require.Equal(t, want.IsMontgomery, got.IsMontgomery)
	for component := range got.Value {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, want.Value[component].Coeffs[limb], got.Value[component].Coeffs[limb])
		}
	}
}

func TestFastTraceGapOneInPlaceAndOutOfPlace(t *testing.T) {
	params := traceTestParameters(t)
	eval := NewEvaluator(params)
	in := newTraceCiphertext(params, 4, false, true)
	want := traceReference(params, in, 3)
	out := newTraceCiphertext(params, 4, false, true)
	metadata := *in.MetaData
	require.NoError(t, eval.Trace(in, 3, out))
	requireTraceMatches(t, want, out)
	require.Equal(t, metadata, *out.MetaData)

	alias := in.CopyNew()
	require.NoError(t, eval.Trace(alias, 3, alias))
	requireTraceMatches(t, want, alias)
}

func TestFastTraceNontrivialAndLogNZero(t *testing.T) {
	params := traceTestParameters(t)
	eval := NewEvaluator(params)
	for _, logN := range []int{1, 0} {
		in := newTraceCiphertext(params, 4, false, false)
		want := traceReference(params, in, logN)
		out := newTraceCiphertext(params, 4, false, true)
		require.NoError(t, eval.Trace(in, logN, out))
		requireTraceMatches(t, want, out)
	}
	require.Len(t, eval.automorphismIndexCache, 4)
}

func TestFastTraceMontgomeryAndPoisonedResidues(t *testing.T) {
	params := traceTestParameters(t)
	eval := NewEvaluator(params)
	for _, montgomery := range []bool{false, true} {
		in := newTraceCiphertext(params, 4, montgomery, true)
		want := traceReference(params, in, 1)
		dormant := cloneDormant(in, 2)
		require.NoError(t, eval.Trace(in, 1, in))
		requireTraceMatches(t, want, in)
		require.Equal(t, dormant, cloneDormant(in, 2))
	}
}

func TestFastTraceValidation(t *testing.T) {
	params := traceTestParameters(t)
	eval := NewEvaluator(params)
	in := newTraceCiphertext(params, 4, false, false)
	_, _ = eval, in
	require.Error(t, eval.Trace(nil, 1, in))
	require.Error(t, eval.Trace(in, 1, nil))
	require.Error(t, eval.Trace(in, -1, in))
	require.Error(t, eval.Trace(in, params.LogN(), in))

	low := newTraceCiphertext(params, 0, false, false)
	require.Error(t, eval.Trace(low, 1, low))
	nonNTT := newTraceCiphertext(params, 4, false, false)
	nonNTT.IsNTT = false
	require.Error(t, eval.Trace(nonNTT, 1, nonNTT))

	mismatchLevel := newTraceCiphertext(params, 3, false, false)
	require.Error(t, eval.Trace(in, 1, mismatchLevel))
	mismatchRepresentation := newTraceCiphertext(params, 4, true, false)
	require.Error(t, eval.Trace(in, 1, mismatchRepresentation))

	wrongDegree := ckks.NewCiphertext(params, 2, 4)
	wrongDegree.IsNTT = true
	require.Error(t, eval.Trace(wrongDegree, 1, wrongDegree))

	ci, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 4, LogQ: []int{50, 50, 50}, RingType: ring.ConjugateInvariant})
	require.NoError(t, err)
	ciEval := NewEvaluator(ci)
	require.Error(t, ciEval.Trace(ckks.NewCiphertext(ci, 1, 1), 1, ckks.NewCiphertext(ci, 1, 1)))
}

func traceBenchmarkParameters(b *testing.B, logN int) ckks.Parameters {
	logQ := []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56}
	if logN == 16 {
		logQ[0], logQ[1] = 54, 38
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: logN, LogQ: logQ, LogDefaultScale: 30})
	if err != nil {
		b.Fatal(err)
	}
	return params
}

func BenchmarkFastTraceLogN13(b *testing.B) {
	benchmarkFastTrace(b, 13)
}

func BenchmarkFastTraceLogN16(b *testing.B) {
	benchmarkFastTrace(b, 16)
}

func benchmarkFastTrace(b *testing.B, logN int) {
	params := traceBenchmarkParameters(b, logN)
	eval := NewEvaluator(params)
	in := newTraceCiphertext(params, params.MaxLevel(), false, true)
	out := newTraceCiphertext(params, params.MaxLevel(), false, true)
	if err := eval.Trace(in, 1, out); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if err := eval.Trace(in, 1, out); err != nil {
			b.Fatal(err)
		}
	}
}
