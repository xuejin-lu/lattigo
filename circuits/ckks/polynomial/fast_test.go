package polynomial

import (
	"math/bits"
	"testing"

	"github.com/stretchr/testify/require"

	commonpolynomial "github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func fastPolynomialTestParameters(t *testing.T) ckks.Parameters {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 35, 35, 35, 35, 35, 35, 35},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	return params
}

func fastPolynomialTestCiphertext(params ckks.Parameters, offset uint64) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, params.MaxLevel())
	ct.IsNTT = true
	ct.IsMontgomery = true
	ct.Scale = params.DefaultScale()
	for d := range ct.Value {
		for limb, subring := range params.RingQ().SubRings[:ct.Level()+1] {
			for i := range ct.Value[d].Coeffs[limb] {
				ct.Value[d].Coeffs[limb][i] = (offset + uint64(17*d+13*limb+i)) % subring.Modulus
			}
		}
	}
	return ct
}

func fastPolynomialTestPoly(degree int) bignum.Polynomial {
	coeffs := make([]*bignum.Complex, degree+1)
	for i := range coeffs {
		coeffs[i] = bignum.ToComplex(int64(i+1), 128)
	}
	poly := bignum.NewPolynomial(bignum.Chebyshev, coeffs, nil)
	poly.IsEven = false
	poly.IsOdd = false
	return poly
}

func TestFastPolynomialEvaluatorSmokeAndReuse(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	poly := fastPolynomialTestPoly(3)

	input := fastPolynomialTestCiphertext(params, 7)
	got, err := eval.Evaluate(input, poly, params.DefaultScale())
	require.NoError(t, err)
	require.Equal(t, 1, got.Degree())
	require.Equal(t, input.IsNTT, got.IsNTT)
	require.Equal(t, input.IsMontgomery, got.IsMontgomery)
	require.GreaterOrEqual(t, got.Level(), 1)

	first := got.CopyNew()
	secondInput := fastPolynomialTestCiphertext(params, 101)
	second, err := eval.Evaluate(secondInput, poly, params.DefaultScale())
	require.NoError(t, err)
	require.NotEqual(t, first.Value[0].Coeffs[0], second.Value[0].Coeffs[0])
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, first.Value[d].Coeffs[limb], got.Value[d].Coeffs[limb])
		}
	}
	require.NotSame(t, first, second)
}

func TestFastPolynomialPowerBasisQ01Oracle(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 7)
	poly := fastPolynomialTestPoly(3)

	_, err := eval.Evaluate(input, poly, params.DefaultScale())
	require.NoError(t, err)

	got := eval.workspace.powers[2]
	reference := ckks.NewCiphertext(params, 1, input.Level())
	reference.IsNTT = input.IsNTT
	reference.IsMontgomery = input.IsMontgomery
	fastEval := eval.Evaluator
	require.NoError(t, fastEval.MulRelin(input, input, reference))
	require.NoError(t, fastEval.Add(reference, reference, reference))
	require.NoError(t, fastEval.Add(reference, -1, reference))
	require.NoError(t, fastEval.Rescale(reference, reference))

	require.Equal(t, reference.Level(), got.Level())
	require.Equal(t, reference.Scale, got.Scale)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, reference.Value[d].Coeffs[limb], got.Value[d].Coeffs[limb])
		}
	}
}

func TestFastPolynomialPlannerMetadata(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 11)
	poly := fastPolynomialTestPoly(7)
	targetScale := params.DefaultScale()

	got, err := eval.Evaluate(input, poly, targetScale)
	require.NoError(t, err)

	commonPoly := commonpolynomial.NewPolynomial(poly)
	levelsConsumed := params.LevelsConsumedPerRescaling()
	sim := simEvaluator{params: params, levelsConsumedPerRescaling: levelsConsumed}
	simPowers := commonpolynomial.SimPowerBasis{1: &commonpolynomial.SimOperand{Level: input.Level(), Scale: input.Scale}}
	logDegree := bits.Len64(uint64(commonPoly.Degree()))
	simPowers.GenPower(eval.Evaluator.GetRLWEParameters(), 1<<logDegree, sim)
	logSplit := bignum.OptimalSplit(logDegree)
	for i := (1 << logSplit) - 1; i > 2; i-- {
		simPowers.GenPower(eval.Evaluator.GetRLWEParameters(), i, sim)
	}
	plan := commonPoly.PatersonStockmeyerPolynomial(eval.Evaluator.GetRLWEParameters(), input.Level(), input.Scale, targetScale, sim)

	type meta struct {
		degree int
		level  int
		scale  rlwe.Scale
	}
	steps := make([]*meta, len(plan.Value))
	for i := range plan.Value {
		steps[len(plan.Value)-i-1] = &meta{degree: plan.Value[i].Degree(), level: plan.Value[i].Level, scale: plan.Value[i].Scale}
	}
	for len(steps) != 1 {
		giant := make([]int, len(steps))
		for i := 0; i < len(steps); i++ {
			if i == len(steps)-1 {
				giant[i] = 2
			} else if steps[i].degree == steps[i+1].degree {
				giant[i] = 1
				i++
			}
		}
		for i := 0; i < len(steps); i++ {
			if giant[i] == 2 {
				steps[i].degree = steps[i-1].degree
			} else if giant[i] == 1 {
				left, right := steps[i], steps[i+1]
				degree := 1 << bits.Len64(uint64(left.degree))
				rescaled := &commonpolynomial.SimOperand{Level: right.level, Scale: right.scale}
				sim.Rescale(rescaled)
				xpow := simPowers[degree]
				level := rescaled.Level
				if xpow.Level < level {
					level = xpow.Level
				}
				scale := rescaled.Scale.Mul(xpow.Scale)
				if left.scale.Cmp(scale) > 0 {
					scale = left.scale
				}
				steps[i+1] = &meta{degree: 2*degree - 1, level: level, scale: scale}
				steps[i] = nil
				i++
			}
		}
		compact := steps[:0]
		for _, step := range steps {
			if step != nil {
				compact = append(compact, step)
			}
		}
		steps = compact
	}
	final := &commonpolynomial.SimOperand{Level: steps[0].level, Scale: steps[0].scale}
	sim.Rescale(final)
	require.Equal(t, final.Level, got.Level())
	require.True(t, final.Scale.Equal(got.Scale), "planned scale=%v got=%v", final.Scale.Float64(), got.Scale.Float64())
}

func TestFastPolynomialIgnoresDormantResidues(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	poly := fastPolynomialTestPoly(3)
	inputA := fastPolynomialTestCiphertext(params, 19)
	inputB := inputA.CopyNew()
	for d := range inputB.Value {
		for limb := 2; limb <= inputB.Level(); limb++ {
			for i := range inputB.Value[d].Coeffs[limb] {
				inputB.Value[d].Coeffs[limb][i] ^= uint64(0x5a5a5a5a) + uint64(31*d+limb+i)
			}
		}
	}

	gotA, err := NewFastEvaluator(params, nil).Evaluate(inputA, poly, params.DefaultScale())
	require.NoError(t, err)
	gotB, err := NewFastEvaluator(params, nil).Evaluate(inputB, poly, params.DefaultScale())
	require.NoError(t, err)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, gotA.Value[d].Coeffs[limb], gotB.Value[d].Coeffs[limb])
		}
	}
}

func TestFastPolynomialRepresentativeDegree30(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            6,
		LogQ:            []int{55, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 23)
	got, err := eval.Evaluate(input, fastPolynomialTestPoly(30), params.DefaultScale())
	require.NoError(t, err)
	require.Equal(t, 1, got.Degree())
	require.Equal(t, input.IsNTT, got.IsNTT)
	require.Equal(t, input.IsMontgomery, got.IsMontgomery)
	require.GreaterOrEqual(t, got.Level(), 1)
}

func TestFastPolynomialEvaluatorValidation(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 1)
	poly := fastPolynomialTestPoly(3)
	_, err := eval.Evaluate(input, bignum.Polynomial{}, params.DefaultScale())
	require.Error(t, err)

	_, err = eval.Evaluate(nil, poly, params.DefaultScale())
	require.Error(t, err)

	nonNTT := input.CopyNew()
	nonNTT.IsNTT = false
	_, err = eval.Evaluate(nonNTT, poly, params.DefaultScale())
	require.Error(t, err)

	nonMontgomery := input.CopyNew()
	nonMontgomery.IsMontgomery = false
	_, err = eval.Evaluate(nonMontgomery, poly, params.DefaultScale())
	require.Error(t, err)

	levelZero := input.CopyNew()
	levelZero.Resize(1, 0)
	_, err = eval.Evaluate(levelZero, poly, params.DefaultScale())
	require.Error(t, err)

	degreeTwo := input.CopyNew()
	degreeTwo.Resize(2, degreeTwo.Level())
	_, err = eval.Evaluate(degreeTwo, poly, params.DefaultScale())
	require.Error(t, err)

	monomial := poly.Clone()
	monomial.Basis = bignum.Monomial
	_, err = eval.Evaluate(input, monomial, params.DefaultScale())
	require.Error(t, err)

	wrongRing, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 35, 35, 35},
		LogDefaultScale: 30,
		RingType:        ring.ConjugateInvariant,
	})
	require.NoError(t, err)
	wrongRingInput := fastPolynomialTestCiphertext(wrongRing, 1)
	_, err = NewFastEvaluator(wrongRing, nil).Evaluate(wrongRingInput, poly, wrongRing.DefaultScale())
	require.Error(t, err)
}
