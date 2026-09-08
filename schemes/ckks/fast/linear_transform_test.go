package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

func newFastLinearTransform(params ckks.Parameters, level int) lintrans.LinearTransformation {
	lt := lintrans.NewLinearTransformation(params, lintrans.Parameters{
		DiagonalsIndexList:        []int{0, 1, -1},
		LevelQ:                    level,
		LevelP:                    0,
		Scale:                     rlwe.NewScale(2),
		LogDimensions:             ring.Dimensions{Cols: 3},
		LogBabyStepGiantStepRatio: -1,
	})
	lt.IsNTT = true
	lt.IsMontgomery = true

	for diagonal, plaintext := range lt.Vec {
		for limb := 0; limb < 2; limb++ {
			for i := range plaintext.Q.Coeffs[limb] {
				plaintext.Q.Coeffs[limb][i] = uint64(3+diagonal+limb+i) % params.RingQ().SubRings[limb].Modulus
			}
			params.RingQ().SubRings[limb].NTT(plaintext.Q.Coeffs[limb], plaintext.Q.Coeffs[limb])
			params.RingQ().SubRings[limb].MForm(plaintext.Q.Coeffs[limb], plaintext.Q.Coeffs[limb])
		}
		lt.Vec[diagonal] = plaintext
	}
	return lt
}

func referenceFastLinearTransform(params ckks.Parameters, ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation, ctOut *rlwe.Ciphertext) {
	r := params.RingQ().AtLevel(ctIn.Level())
	acc0 := ring.NewPoly(params.N(), 1)
	acc1 := ring.NewPoly(params.N(), 1)
	term0 := ring.NewPoly(params.N(), 1)
	term1 := ring.NewPoly(params.N(), 1)
	rot0 := ring.NewPoly(params.N(), 1)
	rot1 := ring.NewPoly(params.N(), 1)
	first := true

	for diagonal, plaintext := range matrix.Vec {
		diagonal &= (1 << matrix.LogDimensions.Cols) - 1
		galEl := params.GaloisElement(diagonal)
		index, err := ring.AutomorphismNTTIndex(r.N(), r.NthRoot(), galEl)
		if err != nil {
			panic(err)
		}
		in0 := r.NewPoly()
		in1 := r.NewPoly()
		out0 := r.NewPoly()
		out1 := r.NewPoly()
		copyQ01(ctIn.Value[0], in0)
		copyQ01(ctIn.Value[1], in1)
		r.AutomorphismNTTWithIndex(in0, index, out0)
		r.AutomorphismNTTWithIndex(in1, index, out1)
		copyQ01(out0, rot0)
		copyQ01(out1, rot1)
		for limb := 0; limb < 2; limb++ {
			r.SubRings[limb].MulCoeffsMontgomery(plaintext.Q.Coeffs[limb], rot0.Coeffs[limb], term0.Coeffs[limb])
			r.SubRings[limb].MulCoeffsMontgomery(plaintext.Q.Coeffs[limb], rot1.Coeffs[limb], term1.Coeffs[limb])
		}
		if first {
			copyQ01(term0, acc0)
			copyQ01(term1, acc1)
			first = false
		} else {
			for limb := 0; limb < 2; limb++ {
				r.SubRings[limb].Add(term0.Coeffs[limb], acc0.Coeffs[limb], acc0.Coeffs[limb])
				r.SubRings[limb].Add(term1.Coeffs[limb], acc1.Coeffs[limb], acc1.Coeffs[limb])
			}
		}
	}
	copyQ01(acc0, ctOut.Value[0])
	copyQ01(acc1, ctOut.Value[1])
}

func TestFastLinearTransformMatchesRingReference(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)

	for _, level := range []int{1, 2, 3} {
		matrix := newFastLinearTransform(params, level)
		in := ckks.NewCiphertext(params, 1, level)
		out := ckks.NewCiphertext(params, 1, level)
		want := ckks.NewCiphertext(params, 1, level)
		fillFastCKKSCiphertext(in, params, nil, 19)
		fillFastCKKSCiphertext(out, params, nil, 101)
		in.IsNTT, in.IsMontgomery = true, true
		out.IsNTT, out.IsMontgomery = true, true
		in.Scale = rlwe.NewScale(7)
		in.IsBatched = true
		in.IsBitReversed = true
		in.LogDimensions = ring.Dimensions{Rows: 1, Cols: 3}
		*out.MetaData = *in.MetaData
		want.IsNTT, want.IsMontgomery = true, true
		r := params.RingQ().AtLevel(level)
		for d := 0; d <= 1; d++ {
			for limb := 0; limb < 2; limb++ {
				r.SubRings[limb].NTT(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
				r.SubRings[limb].MForm(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
			}
		}
		outDormant := cloneDormant(out, 2)
		want.MetaData = in.MetaData.CopyNew()
		want.MetaData.Scale = in.Scale.Mul(matrix.Scale)
		referenceFastLinearTransform(params, in, matrix, want)

		require.NoError(t, eval.LinearTransform(in, matrix, out))
		require.Equal(t, want.Value[0].Coeffs[:2], out.Value[0].Coeffs[:2])
		require.Equal(t, want.Value[1].Coeffs[:2], out.Value[1].Coeffs[:2])
		require.Equal(t, outDormant[:2], cloneDormant(out, 2))
		require.Equal(t, in.Scale.Mul(matrix.Scale), out.Scale)
		metadata := *in.MetaData
		metadata.Scale = out.Scale
		require.Equal(t, metadata, *out.MetaData)
	}
}

func newFastBSGSLinearTransform(params ckks.Parameters, level int) lintrans.LinearTransformation {
	lt := lintrans.NewLinearTransformation(params, lintrans.Parameters{
		DiagonalsIndexList:        []int{0, 1, 2, 3, 4, 5, 6},
		LevelQ:                    level,
		LevelP:                    0,
		Scale:                     rlwe.NewScale(3),
		LogDimensions:             ring.Dimensions{Cols: 3},
		LogBabyStepGiantStepRatio: 1,
	})
	lt.IsNTT = true
	lt.IsMontgomery = true
	for diagonal, plaintext := range lt.Vec {
		for limb := 0; limb < 2; limb++ {
			for i := range plaintext.Q.Coeffs[limb] {
				plaintext.Q.Coeffs[limb][i] = uint64(11+diagonal+3*limb+i) % params.RingQ().SubRings[limb].Modulus
			}
			params.RingQ().SubRings[limb].NTT(plaintext.Q.Coeffs[limb], plaintext.Q.Coeffs[limb])
			params.RingQ().SubRings[limb].MForm(plaintext.Q.Coeffs[limb], plaintext.Q.Coeffs[limb])
		}
		lt.Vec[diagonal] = plaintext
	}
	return lt
}

func referenceFastBSGS(params ckks.Parameters, ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation, ctOut *rlwe.Ciphertext) {
	ringQ := params.RingQ().AtLevel(ctIn.Level())
	index, _, _ := matrix.BSGSIndex()
	acc0, acc1 := ring.NewPoly(params.N(), 1), ring.NewPoly(params.N(), 1)
	inner0, inner1 := ring.NewPoly(params.N(), 1), ring.NewPoly(params.N(), 1)
	rot0, rot1 := ring.NewPoly(params.N(), 1), ring.NewPoly(params.N(), 1)
	term0, term1 := ring.NewPoly(params.N(), 1), ring.NewPoly(params.N(), 1)
	outer0, outer1 := ring.NewPoly(params.N(), 1), ring.NewPoly(params.N(), 1)
	firstOuter := true
	slots := 1 << matrix.LogDimensions.Cols
	for _, j := range utils.GetSortedKeys(index) {
		firstInner := true
		for _, i := range index[j] {
			if i == 0 {
				copyQ01(ctIn.Value[0], rot0)
				copyQ01(ctIn.Value[1], rot1)
			} else {
				if err := FastAutomorphism(ringQ, ctIn.Value[0], rot0, params.GaloisElement(i), true); err != nil {
					panic(err)
				}
				if err := FastAutomorphism(ringQ, ctIn.Value[1], rot1, params.GaloisElement(i), true); err != nil {
					panic(err)
				}
			}
			pt := matrix.Vec[j+i]
			if pt.Q.N() == 0 {
				pt = matrix.Vec[j+i-slots]
			}
			fastPlaintextMul(ringQ, pt.Q, rot0, term0)
			fastPlaintextMul(ringQ, pt.Q, rot1, term1)
			if firstInner {
				copyQ01(term0, inner0)
				copyQ01(term1, inner1)
				firstInner = false
			} else {
				addQ01(ringQ, term0, inner0)
				addQ01(ringQ, term1, inner1)
			}
		}
		if j == 0 {
			copyQ01(inner0, outer0)
			copyQ01(inner1, outer1)
		} else {
			if err := FastAutomorphism(ringQ, inner0, outer0, params.GaloisElement(j), true); err != nil {
				panic(err)
			}
			if err := FastAutomorphism(ringQ, inner1, outer1, params.GaloisElement(j), true); err != nil {
				panic(err)
			}
		}
		if firstOuter {
			copyQ01(outer0, acc0)
			copyQ01(outer1, acc1)
			firstOuter = false
		} else {
			addQ01(ringQ, outer0, acc0)
			addQ01(ringQ, outer1, acc1)
		}
	}
	copyQ01(acc0, ctOut.Value[0])
	copyQ01(acc1, ctOut.Value[1])
}

func TestFastLinearTransformBSGSAndLowerLevel(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	matrix := newFastBSGSLinearTransform(params, 3)
	require.NotZero(t, matrix.N1)
	in := ckks.NewCiphertext(params, 1, 2)
	out := ckks.NewCiphertext(params, 1, 2)
	in.IsNTT, in.IsMontgomery, in.IsBatched = true, true, true
	fillFastCKKSCiphertext(in, params, nil, 137)
	r := params.RingQ().AtLevel(in.Level())
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			r.SubRings[limb].NTT(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
			r.SubRings[limb].MForm(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
		}
	}
	out.IsNTT, out.IsMontgomery, out.IsBatched = true, true, true
	want := ckks.NewCiphertext(params, 1, 2)
	want.IsNTT, want.IsMontgomery, want.IsBatched = true, true, true
	referenceFastBSGS(params, in, matrix, want)
	require.NoError(t, eval.LinearTransform(in, matrix, out))
	require.Equal(t, want.Value[0].Coeffs[:2], out.Value[0].Coeffs[:2])
	require.Equal(t, want.Value[1].Coeffs[:2], out.Value[1].Coeffs[:2])
	require.Equal(t, in.Scale.Mul(matrix.Scale), out.Scale)
}

func TestFastLinearTransformBSGSPrecomputesBabyRotationsOnce(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	matrix := newFastBSGSLinearTransform(params, 3)
	in := ckks.NewCiphertext(params, 1, 3)
	out := ckks.NewCiphertext(params, 1, 3)
	in.IsNTT, in.IsMontgomery, in.IsBatched = true, true, true
	out.IsNTT, out.IsMontgomery, out.IsBatched = true, true, true
	fillFastCKKSCiphertext(in, params, nil, 173)
	r := params.RingQ().AtLevel(in.Level())
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			r.SubRings[limb].NTT(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
			r.SubRings[limb].MForm(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
		}
	}

	index, _, rotN2 := matrix.BSGSIndex()
	uniqueNonZero := 0
	for _, i := range rotN2 {
		if i != 0 {
			uniqueNonZero++
		}
	}
	repeatedNonZero := 0
	for _, rotations := range index {
		for _, i := range rotations {
			if i != 0 {
				repeatedNonZero++
			}
		}
	}

	require.NoError(t, eval.LinearTransform(in, matrix, out))
	t.Logf("BSGS baby rotations: %d unique nonzero versus %d repeated nonzero uses", uniqueNonZero, repeatedNonZero)
	require.Equal(t, uniqueNonZero, eval.lastBSGSBabyRotations)
	require.Greater(t, repeatedNonZero, uniqueNonZero)
}

func TestFastLinearTransformIgnoresDormantAndSupportsAlias(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	level := 3
	matrix := newFastLinearTransform(params, level)
	in := ckks.NewCiphertext(params, 1, level)
	in.IsNTT, in.IsMontgomery = true, true
	fillFastCKKSCiphertext(in, params, nil, 29)
	r := params.RingQ().AtLevel(level)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			r.SubRings[limb].NTT(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
			r.SubRings[limb].MForm(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
		}
	}
	baseline := in.CopyNew()
	garbage := in.CopyNew()
	for d := 0; d <= 1; d++ {
		for limb := 2; limb <= level; limb++ {
			for i := range garbage.Value[d].Coeffs[limb] {
				garbage.Value[d].Coeffs[limb][i] = ^uint64(0) - uint64(i+limb+d)
			}
		}
	}
	matrixGarbage := matrix
	matrixGarbage.Vec = make(map[int]ringqp.Poly, len(matrix.Vec))
	for diagonal, plaintext := range matrix.Vec {
		garbagePlaintext := *plaintext.CopyNew()
		for limb := 2; limb < len(garbagePlaintext.Q.Coeffs); limb++ {
			for i := range garbagePlaintext.Q.Coeffs[limb] {
				garbagePlaintext.Q.Coeffs[limb][i] = ^uint64(0) - uint64(i+limb)
			}
		}
		for limb := range garbagePlaintext.P.Coeffs {
			for i := range garbagePlaintext.P.Coeffs[limb] {
				garbagePlaintext.P.Coeffs[limb][i] = ^uint64(0) - uint64(i+limb+17)
			}
		}
		matrixGarbage.Vec[diagonal] = garbagePlaintext
	}

	outA := ckks.NewCiphertext(params, 1, level)
	outB := ckks.NewCiphertext(params, 1, level)
	fillFastCKKSCiphertext(outA, params, nil, 71)
	fillFastCKKSCiphertext(outB, params, nil, 73)
	dormantA := cloneDormant(outA, 2)
	dormantB := cloneDormant(outB, 2)
	require.NoError(t, eval.LinearTransform(baseline, matrix, outA))
	require.NoError(t, eval.LinearTransform(garbage, matrixGarbage, outB))
	require.Equal(t, outA.Value[0].Coeffs[:2], outB.Value[0].Coeffs[:2])
	require.Equal(t, outA.Value[1].Coeffs[:2], outB.Value[1].Coeffs[:2])
	require.Equal(t, dormantA[:2], cloneDormant(outA, 2))
	require.Equal(t, dormantB[:2], cloneDormant(outB, 2))

	alias := baseline.CopyNew()
	want := outA.CopyNew()
	require.NoError(t, eval.LinearTransform(baseline, matrix, want))
	require.NoError(t, eval.LinearTransform(alias, matrix, alias))
	require.Equal(t, want.Value[0].Coeffs[:2], alias.Value[0].Coeffs[:2])
	require.Equal(t, want.Value[1].Coeffs[:2], alias.Value[1].Coeffs[:2])
}

func TestFastLinearTransformRejectsLevelMismatch(t *testing.T) {
	params := testFastCKKSParameters(t)
	matrix := newFastLinearTransform(params, 2)
	ct := ckks.NewCiphertext(params, 1, 2)
	ct.IsNTT, ct.IsMontgomery = true, true
	eval := NewEvaluator(params)
	matrix = newFastLinearTransform(params, 1)
	require.Error(t, eval.LinearTransform(ct, matrix, ct.CopyNew()))

	coefficientMatrix := newFastLinearTransform(params, 2)
	coefficientMatrix.IsNTT, coefficientMatrix.IsMontgomery = false, false
	coefficientInput := ckks.NewCiphertext(params, 1, 2)
	coefficientInput.IsNTT, coefficientInput.IsMontgomery = false, false
	require.Error(t, eval.LinearTransform(coefficientInput, coefficientMatrix, coefficientInput.CopyNew()))
}

func BenchmarkFastLinearTransform(b *testing.B) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{50, 50, 50, 50, 50, 50, 50, 50, 50},
		LogP:            []int{50},
		LogDefaultScale: 40,
	})
	if err != nil {
		b.Fatal(err)
	}
	level := params.MaxLevel()
	matrix := newFastLinearTransform(params, level)
	in := ckks.NewCiphertext(params, 1, level)
	out := ckks.NewCiphertext(params, 1, level)
	in.IsNTT, in.IsMontgomery = true, true
	out.IsNTT, out.IsMontgomery = true, true
	fillFastCKKSCiphertextBenchmark(in, params, 17)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			params.RingQ().SubRings[limb].NTT(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
			params.RingQ().SubRings[limb].MForm(in.Value[d].Coeffs[limb], in.Value[d].Coeffs[limb])
		}
	}
	eval := NewEvaluator(params)
	kgen := rlwe.NewKeyGenerator(params.Parameters)
	sk := kgen.GenSecretKeyNew()
	galEls := matrix.GaloisElements(params)
	gks := kgen.GenGaloisKeysNew(galEls, sk)
	standardEval := ckks.NewEvaluator(params, rlwe.NewMemEvaluationKeySet(nil, gks...))
	standardLTEval := lintrans.NewEvaluator(standardEval)
	standardIn := ckks.NewCiphertext(params, 1, level)
	standardIn.IsNTT = false
	fillFastCKKSCiphertextBenchmark(standardIn, params, 17)
	standardOut := ckks.NewCiphertext(params, 1, level)
	for d := 0; d <= 1; d++ {
		params.RingQ().NTT(standardIn.Value[d], standardIn.Value[d])
	}
	standardIn.IsNTT, standardIn.IsMontgomery = true, false
	standardOut.IsNTT, standardOut.IsMontgomery = true, true
	b.ReportMetric(float64(len(matrix.Vec)), "diagonals")
	b.Run("FastQ01", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err := eval.LinearTransform(in, matrix, out); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("StandardFullRNS", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if err := standardLTEval.EvaluateMany(standardIn, []lintrans.LinearTransformation{matrix}, []*rlwe.Ciphertext{standardOut}); err != nil {
				b.Fatal(err)
			}
		}
	})
}
