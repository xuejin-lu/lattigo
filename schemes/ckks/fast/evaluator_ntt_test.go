package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func newNTTFastCiphertext(params ckks.Parameters, degree, level int, offset uint64) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, degree, level)
	fillFastCKKSCiphertext(ct, params, nil, offset)
	for d := range ct.Value {
		params.RingQ().AtLevel(level).NTT(ct.Value[d], ct.Value[d])
	}
	ct.IsNTT = true
	return ct
}

func poisonDormant(ct *rlwe.Ciphertext, salt uint64) {
	for d := range ct.Value {
		for limb := 2; limb <= ct.Level(); limb++ {
			for i := range ct.Value[d].Coeffs[limb] {
				ct.Value[d].Coeffs[limb][i] = ^uint64(0) - salt - uint64(31*d+7*limb+i)
			}
		}
	}
}

func requireQ01Equal(t *testing.T, want, got *rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, want.Degree(), got.Degree())
	for d := range want.Value {
		require.Equal(t, want.Value[d].Coeffs[:2], got.Value[d].Coeffs[:2], "component %d", d)
	}
}

func TestFastEvaluatorAddSubNTT(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	for _, tc := range []struct {
		name string
		sub  bool
	}{{"Add", false}, {"Sub", true}} {
		t.Run(tc.name, func(t *testing.T) {
			a := newNTTFastCiphertext(params, 2, 3, 7)
			b := newNTTFastCiphertext(params, 1, 2, 19)
			out := ckks.NewCiphertext(params, 2, 2)
			out.IsNTT = true
			poisonDormant(out, 101)
			dormant := cloneDormant(out, 2)
			want := ckks.NewCiphertext(params, 2, 2)
			standard := ckks.NewEvaluator(params, nil)
			if tc.sub {
				require.NoError(t, standard.Sub(a, b, want))
				require.NoError(t, eval.Sub(a, b, out))
			} else {
				require.NoError(t, standard.Add(a, b, want))
				require.NoError(t, eval.Add(a, b, out))
			}
			requireQ01Equal(t, want, out)
			require.Equal(t, 2, out.Level())
			require.Equal(t, a.Scale, out.Scale)
			require.Equal(t, dormant, cloneDormant(out, 2))
			aPoison, bPoison := a.CopyNew(), b.CopyNew()
			poisonDormant(aPoison, 211)
			poisonDormant(bPoison, 307)
			poisonOut := ckks.NewCiphertext(params, 2, 2)
			poisonOut.IsNTT = true
			if tc.sub {
				require.NoError(t, eval.Sub(aPoison, bPoison, poisonOut))
			} else {
				require.NoError(t, eval.Add(aPoison, bPoison, poisonOut))
			}
			requireQ01Equal(t, out, poisonOut)

			inPlace := a.CopyNew()
			if tc.sub {
				require.NoError(t, eval.Sub(inPlace, b, inPlace))
			} else {
				require.NoError(t, eval.Add(inPlace, b, inPlace))
			}
			requireQ01Equal(t, want, inPlace)
			if tc.sub {
				lowDegree := newNTTFastCiphertext(params, 1, 2, 61)
				highDegree := newNTTFastCiphertext(params, 2, 2, 73)
				wantAlias := ckks.NewCiphertext(params, 2, 2)
				require.NoError(t, standard.Sub(lowDegree, highDegree, wantAlias))
				require.NoError(t, eval.Sub(lowDegree, highDegree, highDegree))
				requireQ01Equal(t, wantAlias, highDegree)
			}

			oneA := newNTTFastCiphertext(params, 1, 3, 37)
			oneB := newNTTFastCiphertext(params, 1, 3, 43)
			var oneOut *rlwe.Ciphertext
			var err error
			if tc.sub {
				oneOut, err = eval.SubNew(oneA, oneB)
			} else {
				oneOut, err = eval.AddNew(oneA, oneB)
			}
			require.NoError(t, err)
			require.Equal(t, 1, oneOut.Degree())
		})
	}
}

func TestFastEvaluatorMulAndMulRelinNTT(t *testing.T) {
	params := testFastCKKSParameters(t)
	standard := ckks.NewEvaluator(params, nil)
	for _, montgomery := range []bool{false, true} {
		t.Run(map[bool]string{false: "NTT", true: "NTTMontgomery"}[montgomery], func(t *testing.T) {
			eval := NewEvaluator(params)
			a := newNTTFastCiphertext(params, 1, 3, 11)
			b := newNTTFastCiphertext(params, 1, 2, 29)
			want := ckks.NewCiphertext(params, 2, 2)
			require.NoError(t, standard.Mul(a, b, want))
			if montgomery {
				for _, ct := range []*rlwe.Ciphertext{a, b} {
					for d := range ct.Value {
						params.RingQ().AtLevel(ct.Level()).MForm(ct.Value[d], ct.Value[d])
					}
					ct.IsMontgomery = true
				}
			}
			out := ckks.NewCiphertext(params, 2, 2)
			out.IsNTT, out.IsMontgomery = true, montgomery
			poisonDormant(out, 53)
			dormant := cloneDormant(out, 2)
			require.NoError(t, eval.Mul(a, b, out))
			require.Equal(t, 2, out.Degree())
			require.Equal(t, a.Scale.Mul(b.Scale), out.Scale)
			if montgomery {
				for d := range out.Value {
					for limb := 0; limb < 2; limb++ {
						params.RingQ().SubRings[limb].IMForm(out.Value[d].Coeffs[limb], out.Value[d].Coeffs[limb])
					}
				}
			}
			requireQ01Equal(t, want, out)
			require.Equal(t, dormant, cloneDormant(out, 2))

			relin := ckks.NewCiphertext(params, 1, 2)
			product := ckks.NewCiphertext(params, 2, 2)
			mulThenRelin := ckks.NewCiphertext(params, 1, 2)
			for _, ct := range []*rlwe.Ciphertext{relin, product, mulThenRelin} {
				ct.IsNTT, ct.IsMontgomery = true, montgomery
			}
			require.NoError(t, eval.MulRelin(a, b, relin))
			require.NoError(t, eval.Mul(a, b, product))
			require.NoError(t, eval.Relinearize(product, mulThenRelin))
			requireQ01Equal(t, mulThenRelin, relin)
			require.Equal(t, 1, relin.Degree())
		})
	}
}

func TestFastEvaluatorMulAliasingSquaringAndPoison(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	a := newNTTFastCiphertext(params, 1, 3, 5)
	b := newNTTFastCiphertext(params, 1, 3, 17)
	aPoison, bPoison := a.CopyNew(), b.CopyNew()
	poisonDormant(aPoison, 91)
	poisonDormant(bPoison, 193)
	want, err := eval.MulNew(a, b)
	require.NoError(t, err)
	got, err := eval.MulNew(aPoison, bPoison)
	require.NoError(t, err)
	requireQ01Equal(t, want, got)
	wantRelin, err := eval.MulRelinNew(a, b)
	require.NoError(t, err)
	gotRelin, err := eval.MulRelinNew(aPoison, bPoison)
	require.NoError(t, err)
	requireQ01Equal(t, wantRelin, gotRelin)

	aAlias := a.CopyNew()
	require.NoError(t, eval.Mul(aAlias, b, aAlias))
	requireQ01Equal(t, want, aAlias)
	bAlias := b.CopyNew()
	require.NoError(t, eval.Mul(a, bAlias, bAlias))
	requireQ01Equal(t, want, bAlias)

	square, err := eval.MulNew(a, a)
	require.NoError(t, err)
	standardSquare := ckks.NewCiphertext(params, 2, 3)
	require.NoError(t, ckks.NewEvaluator(params, nil).Mul(a, a, standardSquare))
	requireQ01Equal(t, standardSquare, square)
}

func TestFastEvaluatorRelinearizePoison(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	a := newNTTFastCiphertext(params, 2, 3, 41)
	b := a.CopyNew()
	poisonDormant(a, 7)
	poisonDormant(b, 99)
	outA := ckks.NewCiphertext(params, 1, 3)
	outB := ckks.NewCiphertext(params, 1, 3)
	outA.IsNTT, outB.IsNTT = true, true
	require.NoError(t, eval.Relinearize(a, outA))
	require.NoError(t, eval.Relinearize(b, outB))
	requireQ01Equal(t, outA, outB)
	require.Equal(t, 1, outA.Degree())
	for d := 0; d < 2; d++ {
		require.Equal(t, a.Value[d].Coeffs[:2], outA.Value[d].Coeffs[:2])
	}
}

func TestFastEvaluatorScalarMulThenAdd(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	standard := ckks.NewEvaluator(params, nil)
	x := newNTTFastCiphertext(params, 1, 3, 23)
	coefficient := bignum.ToComplex(complex(1.25, -0.5), params.EncodingPrecision())
	fastOut := newNTTFastCiphertext(params, 1, 3, 71)
	stdOut := fastOut.CopyNew()
	poisonDormant(fastOut, 117)
	poisonDormant(stdOut, 117)
	dormant := cloneDormant(fastOut, 2)
	require.NoError(t, eval.MulThenAdd(x, coefficient, fastOut))
	require.NoError(t, standard.MulThenAdd(x, coefficient, stdOut))
	requireQ01Equal(t, stdOut, fastOut)
	require.Equal(t, stdOut.Scale, fastOut.Scale)
	require.Equal(t, 3, fastOut.Level())
	require.Equal(t, dormant, cloneDormant(fastOut, 2))
	xPoison := x.CopyNew()
	poisonDormant(xPoison, 229)
	poisonOut := newNTTFastCiphertext(params, 1, 3, 71)
	poisonDormant(poisonOut, 313)
	require.NoError(t, eval.MulThenAdd(xPoison, coefficient, poisonOut))
	requireQ01Equal(t, fastOut, poisonOut)

	integer := bignum.ToComplex(complex(2, 0), params.EncodingPrecision())
	alias := x.CopyNew()
	want := x.CopyNew()
	require.NoError(t, standard.MulThenAdd(x, integer, want))
	require.NoError(t, eval.MulThenAdd(alias, integer, alias))
	requireQ01Equal(t, want, alias)

	nonIntegerAlias := x.CopyNew()
	nonIntegerWant := x.CopyNew()
	require.NoError(t, standard.MulThenAdd(x, coefficient, nonIntegerWant))
	require.NoError(t, eval.MulThenAdd(nonIntegerAlias, coefficient, nonIntegerAlias))
	requireQ01Equal(t, nonIntegerWant, nonIntegerAlias)
}

func TestFastEvaluatorScalarAndPlaintextMul(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	standard := ckks.NewEvaluator(params, nil)
	x := newNTTFastCiphertext(params, 1, 3, 47)
	coefficient := bignum.ToComplex(complex(0.75, 0.25), params.EncodingPrecision())
	fastScalar, err := eval.MulNew(x, coefficient)
	require.NoError(t, err)
	stdScalar := ckks.NewCiphertext(params, 1, 3)
	require.NoError(t, standard.Mul(x, coefficient, stdScalar))
	requireQ01Equal(t, stdScalar, fastScalar)
	require.Equal(t, stdScalar.Scale, fastScalar.Scale)
	xPoison := x.CopyNew()
	poisonDormant(xPoison, 401)
	fastScalarPoison, err := eval.MulNew(xPoison, coefficient)
	require.NoError(t, err)
	requireQ01Equal(t, fastScalar, fastScalarPoison)

	ptSource := newNTTFastCiphertext(params, 0, 2, 59)
	pt := ptSource.Plaintext()
	fastPlain, err := eval.MulNew(x, pt)
	require.NoError(t, err)
	stdPlain := ckks.NewCiphertext(params, 1, 2)
	require.NoError(t, standard.Mul(x, pt, stdPlain))
	requireQ01Equal(t, stdPlain, fastPlain)
	require.Equal(t, stdPlain.Scale, fastPlain.Scale)
	ptPoison := pt.CopyNew()
	for limb := 2; limb <= ptPoison.Level(); limb++ {
		for i := range ptPoison.Value.Coeffs[limb] {
			ptPoison.Value.Coeffs[limb][i] = ^uint64(0) - uint64(limb+i)
		}
	}
	fastPlainPoison, err := eval.MulNew(xPoison, ptPoison)
	require.NoError(t, err)
	requireQ01Equal(t, fastPlain, fastPlainPoison)
}

func TestFastEvaluatorScalarMontgomery(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	standard := ckks.NewEvaluator(params, nil)
	x := newNTTFastCiphertext(params, 1, 3, 83)
	coefficient := bignum.ToComplex(complex(1.5, -0.25), params.EncodingPrecision())
	wantAdd, wantMul := ckks.NewCiphertext(params, 1, 3), ckks.NewCiphertext(params, 1, 3)
	require.NoError(t, standard.Add(x, coefficient, wantAdd))
	require.NoError(t, standard.Mul(x, coefficient, wantMul))
	xMont := x.CopyNew()
	for d := range xMont.Value {
		params.RingQ().AtLevel(xMont.Level()).MForm(xMont.Value[d], xMont.Value[d])
	}
	xMont.IsMontgomery = true
	gotAdd, gotMul := ckks.NewCiphertext(params, 1, 3), ckks.NewCiphertext(params, 1, 3)
	gotAdd.IsMontgomery, gotMul.IsMontgomery = true, true
	require.NoError(t, eval.Add(xMont, coefficient, gotAdd))
	require.NoError(t, eval.Mul(xMont, coefficient, gotMul))
	for _, ct := range []*rlwe.Ciphertext{gotAdd, gotMul} {
		for d := range ct.Value {
			for limb := 0; limb < 2; limb++ {
				params.RingQ().SubRings[limb].IMForm(ct.Value[d].Coeffs[limb], ct.Value[d].Coeffs[limb])
			}
		}
	}
	requireQ01Equal(t, wantAdd, gotAdd)
	requireQ01Equal(t, wantMul, gotMul)
}

func TestFastEvaluatorPowerBasisStyleContract(t *testing.T) {
	params := rescaleTestParameters(t)
	eval := NewEvaluator(params)
	x := newNTTFastCiphertext(params, 1, 3, 13)
	x2, err := eval.MulRelinNew(x, x)
	require.NoError(t, err)
	require.Equal(t, 1, x2.Degree())
	require.Equal(t, x.Scale.Mul(x.Scale), x2.Scale)
	require.NoError(t, eval.Rescale(x2, x2))
	require.Equal(t, 2, x2.Level())
	x3, err := eval.MulNew(x2, x)
	require.NoError(t, err)
	require.Equal(t, 2, x3.Degree())
	require.Equal(t, x2.Scale.Mul(x.Scale), x3.Scale)
	require.NoError(t, eval.Relinearize(x3, x3))
	require.Equal(t, 1, x3.Degree())
	require.NoError(t, eval.Rescale(x3, x3))
	require.Equal(t, 1, x3.Level())
}

func TestFastEvaluatorRemainingPoisonPaths(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	for _, sub := range []bool{false, true} {
		a := newNTTFastCiphertext(params, 1, 3, 31)
		b := a.CopyNew()
		poisonDormant(a, 3)
		poisonDormant(b, 79)
		outA := ckks.NewCiphertext(params, 1, 3)
		outB := ckks.NewCiphertext(params, 1, 3)
		outA.IsNTT, outB.IsNTT = true, true
		if sub {
			require.NoError(t, eval.Sub(a, complex(0.5, 0.25), outA))
			require.NoError(t, eval.Sub(b, complex(0.5, 0.25), outB))
		} else {
			require.NoError(t, eval.Add(a, complex(0.5, 0.25), outA))
			require.NoError(t, eval.Add(b, complex(0.5, 0.25), outB))
		}
		requireQ01Equal(t, outA, outB)
	}

	a := newNTTFastCiphertext(params, 1, 3, 7)
	b := newNTTFastCiphertext(params, 1, 3, 17)
	outA := newNTTFastCiphertext(params, 2, 3, 29)
	outA.Scale = a.Scale.Mul(b.Scale)
	outB := outA.CopyNew()
	a2, b2 := a.CopyNew(), b.CopyNew()
	poisonDormant(a2, 103)
	poisonDormant(b2, 211)
	poisonDormant(outB, 307)
	require.NoError(t, eval.MulThenAdd(a, b, outA))
	require.NoError(t, eval.MulThenAdd(a2, b2, outB))
	requireQ01Equal(t, outA, outB)
}

func TestFastEvaluatorRejectsVectorScalar(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	x := newNTTFastCiphertext(params, 1, 3, 1)
	out := ckks.NewCiphertext(params, 1, 3)
	out.IsNTT = true
	require.Error(t, eval.Add(x, []complex128{1}, out))
}
