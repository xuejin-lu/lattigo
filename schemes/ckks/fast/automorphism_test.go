package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func fastAutomorphismTestParameters(t *testing.T, ringType ring.Type) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{50, 50, 50, 50, 50},
		LogDefaultScale: 40,
		RingType:        ringType,
	})
	require.NoError(t, err)
	return params
}

func setFastAutomorphismMetadata(ct *rlwe.Ciphertext) {
	ct.Scale = rlwe.NewScale(123.5)
	ct.IsBatched = true
	ct.IsBitReversed = true
	ct.LogDimensions = ring.Dimensions{Rows: 1, Cols: 3}
	ct.IsMontgomery = true
}

func fastAutomorphismGaloisElements(params ckks.Parameters) []uint64 {
	return []uint64{
		1,
		params.GaloisElementForRotation(1),
		params.GaloisElementForRotation(-1),
		params.GaloisElementForComplexConjugation(),
		params.RingQ().NthRoot() - 1,
	}
}

func TestFastAutomorphismCoefficientDomain(t *testing.T) {
	params := fastAutomorphismTestParameters(t, ring.Standard)
	eval := NewEvaluator(params)

	for _, level := range []int{1, 2, 4} {
		for _, galEl := range fastAutomorphismGaloisElements(params) {
			in := ckks.NewCiphertext(params, 1, level)
			out := ckks.NewCiphertext(params, 1, level)
			fillFastCKKSCiphertext(in, params, nil, 7)
			fillFastCKKSCiphertext(out, params, nil, 101)
			setFastAutomorphismMetadata(in)
			in.IsNTT = false
			out.IsNTT, out.IsMontgomery = false, true
			outDormant := cloneDormant(out, 2)
			metadata := *in.MetaData

			want := in.CopyNew()
			r := params.RingQ().AtLevel(level)
			for d := 0; d <= 1; d++ {
				r.Automorphism(in.Value[d], galEl, want.Value[d])
			}

			require.NoError(t, eval.Automorphism(in, out, galEl))
			require.Equal(t, want.Value[0].Coeffs[:2], out.Value[0].Coeffs[:2])
			require.Equal(t, want.Value[1].Coeffs[:2], out.Value[1].Coeffs[:2])
			require.Equal(t, metadata, *out.MetaData)
			require.True(t, out.IsMontgomery)
			require.Equal(t, outDormant, cloneDormant(out, 2))

			// The same operation must be safe when output aliases input.
			alias := in.CopyNew()
			require.NoError(t, eval.Automorphism(alias, alias, galEl))
			require.Equal(t, want.Value[0].Coeffs[:2], alias.Value[0].Coeffs[:2])
			require.Equal(t, want.Value[1].Coeffs[:2], alias.Value[1].Coeffs[:2])
		}
	}
}

func TestFastAutomorphismNTTDomain(t *testing.T) {
	params := fastAutomorphismTestParameters(t, ring.Standard)
	eval := NewEvaluator(params)

	for _, level := range []int{1, 2, 4} {
		for _, galEl := range fastAutomorphismGaloisElements(params) {
			in := ckks.NewCiphertext(params, 1, level)
			out := ckks.NewCiphertext(params, 1, level)
			fillFastCKKSCiphertext(in, params, nil, 17)
			fillFastCKKSCiphertext(out, params, nil, 201)
			setFastAutomorphismMetadata(in)
			in.IsNTT, out.IsNTT = true, true
			out.IsMontgomery = true
			r := params.RingQ().AtLevel(level)
			for d := 0; d <= 1; d++ {
				for limb := 0; limb < 2; limb++ {
					r.SubRings[limb].NTT(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
				}
			}
			outDormant := cloneDormant(out, 2)

			want := in.CopyNew()
			index, err := ring.AutomorphismNTTIndex(r.N(), r.NthRoot(), galEl)
			require.NoError(t, err)
			for d := 0; d <= 1; d++ {
				r.AutomorphismNTTWithIndex(in.Value[d], index, want.Value[d])
			}

			require.NoError(t, eval.Automorphism(in, out, galEl))
			require.Equal(t, want.Value[0].Coeffs[:2], out.Value[0].Coeffs[:2])
			require.Equal(t, want.Value[1].Coeffs[:2], out.Value[1].Coeffs[:2])
			require.True(t, out.IsNTT)
			require.True(t, out.IsMontgomery)
			require.Equal(t, outDormant, cloneDormant(out, 2))

			alias := in.CopyNew()
			require.NoError(t, eval.Automorphism(alias, alias, galEl))
			require.Equal(t, want.Value[0].Coeffs[:2], alias.Value[0].Coeffs[:2])
			require.Equal(t, want.Value[1].Coeffs[:2], alias.Value[1].Coeffs[:2])
		}
	}
}

func TestFastAutomorphismIgnoresDormantLimbs(t *testing.T) {
	params := fastAutomorphismTestParameters(t, ring.Standard)
	eval := NewEvaluator(params)
	a := ckks.NewCiphertext(params, 1, 4)
	b := ckks.NewCiphertext(params, 1, 4)
	fillFastCKKSCiphertext(a, params, nil, 11)
	fillFastCKKSCiphertext(b, params, nil, 11)
	for d := range b.Value {
		for limb := 2; limb <= b.Level(); limb++ {
			for j := range b.Value[d].Coeffs[limb] {
				b.Value[d].Coeffs[limb][j] = uint64(9000 + d + limb + j)
			}
		}
	}
	a.IsNTT, b.IsNTT = false, false
	outA := ckks.NewCiphertext(params, 1, 4)
	outB := ckks.NewCiphertext(params, 1, 4)
	fillFastCKKSCiphertext(outA, params, nil, 31)
	fillFastCKKSCiphertext(outB, params, nil, 31)
	outA.IsNTT, outB.IsNTT = false, false
	dormantA, dormantB := cloneDormant(outA, 2), cloneDormant(outB, 2)
	galEl := params.GaloisElementForRotation(3)
	require.NoError(t, eval.Automorphism(a, outA, galEl))
	require.NoError(t, eval.Automorphism(b, outB, galEl))
	require.Equal(t, outA.Value[0].Coeffs[:2], outB.Value[0].Coeffs[:2])
	require.Equal(t, outA.Value[1].Coeffs[:2], outB.Value[1].Coeffs[:2])
	require.Equal(t, dormantA, cloneDormant(outA, 2))
	require.Equal(t, dormantB, cloneDormant(outB, 2))
}

func TestFastAutomorphismRejectsConjugateInvariantRing(t *testing.T) {
	params := fastAutomorphismTestParameters(t, ring.ConjugateInvariant)
	eval := NewEvaluator(params)
	in := ckks.NewCiphertext(params, 1, 1)
	out := ckks.NewCiphertext(params, 1, 1)
	require.Error(t, eval.Automorphism(in, out, 1))
	require.Error(t, FastAutomorphism(params.RingQ().AtLevel(1), in.Value[0], out.Value[0], 1, false))
}

func TestFastRotateNew(t *testing.T) {
	params := fastAutomorphismTestParameters(t, ring.Standard)
	eval := NewEvaluator(params)
	in := ckks.NewCiphertext(params, 1, 2)
	fillFastCKKSCiphertext(in, params, nil, 43)
	in.IsNTT = false
	want := ckks.NewCiphertext(params, 1, 2)
	want.IsNTT = false
	require.NoError(t, eval.Automorphism(in, want, params.GaloisElementForRotation(2)))
	got, err := eval.RotateNew(in, 2)
	require.NoError(t, err)
	require.Equal(t, want.Value[0].Coeffs[:2], got.Value[0].Coeffs[:2])
	require.Equal(t, want.Value[1].Coeffs[:2], got.Value[1].Coeffs[:2])
}
