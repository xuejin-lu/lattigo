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
		for limb := range plaintext.Q.Coeffs {
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

func referenceFastLinearTransform(params ckks.Parameters, ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation, ctOut *rlwe.Ciphertext, rows int) {
	r := params.RingQ().AtLevel(rows - 1)
	acc0 := ring.NewPoly(params.N(), rows-1)
	acc1 := ring.NewPoly(params.N(), rows-1)
	term0 := ring.NewPoly(params.N(), rows-1)
	term1 := ring.NewPoly(params.N(), rows-1)
	rot0 := ring.NewPoly(params.N(), rows-1)
	rot1 := ring.NewPoly(params.N(), rows-1)
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
		copyPrefixRowsUnchecked(rows, ctIn.Value[0], in0)
		copyPrefixRowsUnchecked(rows, ctIn.Value[1], in1)
		r.AutomorphismNTTWithIndex(in0, index, out0)
		r.AutomorphismNTTWithIndex(in1, index, out1)
		copyPrefixRowsUnchecked(rows, out0, rot0)
		copyPrefixRowsUnchecked(rows, out1, rot1)
		for row := 0; row < rows; row++ {
			r.SubRings[row].MulCoeffsMontgomery(plaintext.Q.Coeffs[row], rot0.Coeffs[row], term0.Coeffs[row])
			r.SubRings[row].MulCoeffsMontgomery(plaintext.Q.Coeffs[row], rot1.Coeffs[row], term1.Coeffs[row])
		}
		if first {
			copyPrefixRowsUnchecked(rows, term0, acc0)
			copyPrefixRowsUnchecked(rows, term1, acc1)
			first = false
		} else {
			for row := 0; row < rows; row++ {
				r.SubRings[row].Add(term0.Coeffs[row], acc0.Coeffs[row], acc0.Coeffs[row])
				r.SubRings[row].Add(term1.Coeffs[row], acc1.Coeffs[row], acc1.Coeffs[row])
			}
		}
	}
	copyPrefixRowsUnchecked(rows, acc0, ctOut.Value[0])
	copyPrefixRowsUnchecked(rows, acc1, ctOut.Value[1])
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
		referenceFastLinearTransform(params, in, matrix, want, 2)

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

func TestFastLinearTransformQPrefixRowsDirectAndBSGS(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	for _, useBSGS := range []bool{false, true} {
		name := "direct"
		matrix := newFastLinearTransform(params, 3)
		if useBSGS {
			name = "bsgs"
			matrix = newFastBSGSLinearTransform(params, 3)
			require.NotZero(t, matrix.N1)
		} else {
			require.Zero(t, matrix.N1)
		}
		t.Run(name, func(t *testing.T) {
			for _, rows := range []int{3, 4} {
				in := ckks.NewCiphertext(params, 1, 3)
				out := ckks.NewCiphertext(params, 1, 3)
				want := ckks.NewCiphertext(params, 1, 3)
				fillFastCKKSCiphertext(in, params, nil, 211)
				fillFastCKKSCiphertext(out, params, nil, 223)
				in.IsNTT, in.IsMontgomery, in.IsBatched = true, true, true
				out.IsNTT, out.IsMontgomery, out.IsBatched = true, true, true
				for d := 0; d <= 1; d++ {
					for row := 0; row < 4; row++ {
						params.RingQ().SubRings[row].NTT(in.Value[d].Coeffs[row], in.Value[d].Coeffs[row])
						params.RingQ().SubRings[row].MForm(in.Value[d].Coeffs[row], in.Value[d].Coeffs[row])
					}
				}
				untouched := make([][]uint64, 2)
				for d := 0; d <= 1; d++ {
					if rows < len(out.Value[d].Coeffs) {
						untouched[d] = append([]uint64(nil), out.Value[d].Coeffs[rows]...)
					}
				}
				if useBSGS {
					referenceFastBSGS(params, in, matrix, want, rows)
				} else {
					referenceFastLinearTransform(params, in, matrix, want, rows)
				}
				require.NoError(t, eval.LinearTransformQPrefixRows(in, matrix, rows, out))
				for d := 0; d <= 1; d++ {
					require.Equal(t, want.Value[d].Coeffs[:rows], out.Value[d].Coeffs[:rows], "rows=%d component=%d", rows, d)
					if rows < len(out.Value[d].Coeffs) {
						require.Equal(t, untouched[d], out.Value[d].Coeffs[rows], "rows=%d must leave q%d untouched", rows, rows)
					}
				}
				require.Equal(t, in.Scale.Mul(matrix.Scale), out.Scale)
				require.Equal(t, in.IsNTT, out.IsNTT)
				require.Equal(t, in.IsMontgomery, out.IsMontgomery)
				alias := in.CopyNew()
				require.NoError(t, eval.LinearTransformQPrefixRows(alias, matrix, rows, alias))
				for d := 0; d <= 1; d++ {
					require.Equal(t, want.Value[d].Coeffs[:rows], alias.Value[d].Coeffs[:rows], "rows=%d in-place component=%d", rows, d)
				}
			}
		})
	}
}

func TestFastLinearTransformLegacyQ012WidthDirectAndBSGS(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		Q:               []uint64{72057594037616641, 549755731969, 549756026881, 549755486209},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	require.Equal(t, 3, MaintainedLimbCount(&params, 3), "fixture must select legacy q012 authority")
	eval := NewEvaluator(params)
	q3 := params.RingQ().SubRings[3].Modulus

	for _, useBSGS := range []bool{false, true} {
		name := "direct"
		matrix := newFastLinearTransform(params, 3)
		if useBSGS {
			name = "bsgs"
			matrix = newFastBSGSLinearTransform(params, 3)
			require.NotZero(t, matrix.N1)
		} else {
			require.Zero(t, matrix.N1)
		}
		t.Run(name, func(t *testing.T) {
			in := ckks.NewCiphertext(params, 1, 3)
			fillFastCKKSCiphertext(in, params, nil, 307)
			in.IsNTT, in.IsMontgomery, in.IsBatched = true, true, true
			for d := 0; d <= 1; d++ {
				for row := 0; row < 4; row++ {
					params.RingQ().SubRings[row].NTT(in.Value[d].Coeffs[row], in.Value[d].Coeffs[row])
					params.RingQ().SubRings[row].MForm(in.Value[d].Coeffs[row], in.Value[d].Coeffs[row])
				}
			}

			poisonedInput := in.CopyNew()
			for d := 0; d <= 1; d++ {
				for i, residue := range poisonedInput.Value[d].Coeffs[3] {
					poisonedInput.Value[d].Coeffs[3][i] = (residue + uint64(101+i+d)) % q3
				}
			}
			poisonedMatrix := matrix
			poisonedMatrix.Vec = make(map[int]ringqp.Poly, len(matrix.Vec))
			for diagonal, plaintext := range matrix.Vec {
				poisonedPlaintext := *plaintext.CopyNew()
				for i, residue := range poisonedPlaintext.Q.Coeffs[3] {
					poisonedPlaintext.Q.Coeffs[3][i] = (residue + uint64(211+i+diagonal)) % q3
				}
				poisonedMatrix.Vec[diagonal] = poisonedPlaintext
			}

			want := ckks.NewCiphertext(params, 1, 3)
			want.IsNTT, want.IsMontgomery, want.IsBatched = true, true, true
			require.NoError(t, eval.LinearTransformQPrefixRows(in, matrix, 3, want))
			legacyClean := in.CopyNew()
			legacyPoisoned := poisonedInput.CopyNew()
			q3CleanBefore := make([][]uint64, 2)
			q3PoisonedBefore := make([][]uint64, 2)
			for d := 0; d <= 1; d++ {
				q3CleanBefore[d] = append([]uint64(nil), legacyClean.Value[d].Coeffs[3]...)
				q3PoisonedBefore[d] = append([]uint64(nil), legacyPoisoned.Value[d].Coeffs[3]...)
			}
			require.NoError(t, eval.LinearTransform(legacyClean, matrix, legacyClean))
			require.NoError(t, eval.LinearTransform(legacyPoisoned, poisonedMatrix, legacyPoisoned))
			for d := 0; d <= 1; d++ {
				require.Equal(t, want.Value[d].Coeffs[:3], legacyClean.Value[d].Coeffs[:3], "legacy/reference component=%d", d)
				require.Equal(t, legacyClean.Value[d].Coeffs[:3], legacyPoisoned.Value[d].Coeffs[:3], "legacy wrapper consumed q3 poison component=%d", d)
				require.Equal(t, q3CleanBefore[d], legacyClean.Value[d].Coeffs[3], "legacy wrapper must leave q3 untouched component=%d", d)
				require.Equal(t, q3PoisonedBefore[d], legacyPoisoned.Value[d].Coeffs[3], "legacy wrapper must preserve poisoned q3 untouched component=%d", d)
			}
		})
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
		for limb := range plaintext.Q.Coeffs {
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

func referenceFastBSGS(params ckks.Parameters, ctIn *rlwe.Ciphertext, matrix lintrans.LinearTransformation, ctOut *rlwe.Ciphertext, rows int) {
	ringQ := params.RingQ().AtLevel(rows - 1)
	index, _, _ := matrix.BSGSIndex()
	acc0, acc1 := ring.NewPoly(params.N(), rows-1), ring.NewPoly(params.N(), rows-1)
	inner0, inner1 := ring.NewPoly(params.N(), rows-1), ring.NewPoly(params.N(), rows-1)
	rot0, rot1 := ring.NewPoly(params.N(), rows-1), ring.NewPoly(params.N(), rows-1)
	term0, term1 := ring.NewPoly(params.N(), rows-1), ring.NewPoly(params.N(), rows-1)
	outer0, outer1 := ring.NewPoly(params.N(), rows-1), ring.NewPoly(params.N(), rows-1)
	firstOuter := true
	slots := 1 << matrix.LogDimensions.Cols
	for _, j := range utils.GetSortedKeys(index) {
		firstInner := true
		for _, i := range index[j] {
			if i == 0 {
				copyPrefixRowsUnchecked(rows, ctIn.Value[0], rot0)
				copyPrefixRowsUnchecked(rows, ctIn.Value[1], rot1)
			} else {
				referenceFastAutomorphismRows(params, ringQ, ctIn.Value[0], rot0, params.GaloisElement(i), rows)
				referenceFastAutomorphismRows(params, ringQ, ctIn.Value[1], rot1, params.GaloisElement(i), rows)
			}
			pt := matrix.Vec[j+i]
			if pt.Q.N() == 0 {
				pt = matrix.Vec[j+i-slots]
			}
			for row := 0; row < rows; row++ {
				ringQ.SubRings[row].MulCoeffsMontgomery(pt.Q.Coeffs[row], rot0.Coeffs[row], term0.Coeffs[row])
				ringQ.SubRings[row].MulCoeffsMontgomery(pt.Q.Coeffs[row], rot1.Coeffs[row], term1.Coeffs[row])
			}
			if firstInner {
				copyPrefixRowsUnchecked(rows, term0, inner0)
				copyPrefixRowsUnchecked(rows, term1, inner1)
				firstInner = false
			} else {
				for row := 0; row < rows; row++ {
					ringQ.SubRings[row].Add(term0.Coeffs[row], inner0.Coeffs[row], inner0.Coeffs[row])
					ringQ.SubRings[row].Add(term1.Coeffs[row], inner1.Coeffs[row], inner1.Coeffs[row])
				}
			}
		}
		if j == 0 {
			copyPrefixRowsUnchecked(rows, inner0, outer0)
			copyPrefixRowsUnchecked(rows, inner1, outer1)
		} else {
			referenceFastAutomorphismRows(params, ringQ, inner0, outer0, params.GaloisElement(j), rows)
			referenceFastAutomorphismRows(params, ringQ, inner1, outer1, params.GaloisElement(j), rows)
		}
		if firstOuter {
			copyPrefixRowsUnchecked(rows, outer0, acc0)
			copyPrefixRowsUnchecked(rows, outer1, acc1)
			firstOuter = false
		} else {
			for row := 0; row < rows; row++ {
				ringQ.SubRings[row].Add(outer0.Coeffs[row], acc0.Coeffs[row], acc0.Coeffs[row])
				ringQ.SubRings[row].Add(outer1.Coeffs[row], acc1.Coeffs[row], acc1.Coeffs[row])
			}
		}
	}
	copyPrefixRowsUnchecked(rows, acc0, ctOut.Value[0])
	copyPrefixRowsUnchecked(rows, acc1, ctOut.Value[1])
}

func referenceFastAutomorphismRows(params ckks.Parameters, ringQ *ring.Ring, src, dst ring.Poly, galEl uint64, rows int) {
	ringQ = ringQ.AtLevel(rows - 1)
	index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galEl)
	if err != nil {
		panic(err)
	}
	in := ring.NewPoly(params.N(), rows-1)
	copyPrefixRowsUnchecked(rows, src, in)
	result := ring.NewPoly(params.N(), rows-1)
	ringQ.AutomorphismNTTWithIndex(in, index, result)
	copyPrefixRowsUnchecked(rows, result, dst)
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
	referenceFastBSGS(params, in, matrix, want, 2)
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
