package polynomial

import (
	"math"
	"math/big"
	"math/bits"
	"testing"

	"github.com/stretchr/testify/require"

	commonpolynomial "github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/cosine"
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

func fastPolynomialQ012TestParameters(t *testing.T) ckks.Parameters {
	t.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(45, 2*2*16)
	remainingQ, err := generator.NextAlternatingPrimes(5)
	require.NoError(t, err)
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		Q:               append([]uint64{72057594037616641, 549755731969, 549756026881}, remainingQ...),
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

func TestFastPolynomialPowerBasisQPrefixRowsOracle(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 7)
	poly := fastPolynomialTestPoly(3)
	rows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)

	_, err = eval.EvaluateQPrefixRows(input, poly, params.DefaultScale(), rows)
	require.NoError(t, err)

	got := eval.workspace.powers[2]
	reference := ckks.NewCiphertext(params, 1, input.Level())
	reference.IsNTT = input.IsNTT
	reference.IsMontgomery = input.IsMontgomery
	fastEval := eval.Evaluator
	require.NoError(t, fastEval.MulRelinElementQPrefixRows(input, input.El(), rows, reference))
	require.NoError(t, fastEval.AddQPrefixRows(reference, reference, reference, rows))
	require.NoError(t, fastEval.RescaleQPrefixRows(reference, rows, reference))
	rows, err = fastckks.QPrefixWidth(reference.Level())
	require.NoError(t, err)
	require.NoError(t, fastEval.AddScalarQPrefixRows(reference, -1, rows, reference))

	require.Equal(t, reference.Level(), got.Level())
	require.True(t, reference.Scale.Equal(got.Scale), "reference scale=%v got=%v", reference.Scale.Float64(), got.Scale.Float64())
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < rows; limb++ {
			require.Equal(t, reference.Value[d].Coeffs[limb], got.Value[d].Coeffs[limb])
		}
	}
}

func TestFastPolynomialNonPowerOfTwoChebyshevOracle(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            6,
		LogQ:            []int{55, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 13)
	poly := fastPolynomialTestPoly(30)
	rows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)

	_, err = eval.EvaluateQPrefixRows(input, poly, params.DefaultScale(), rows)
	require.NoError(t, err)
	got := eval.workspace.powers[5]
	fastEval := eval.Evaluator

	t2 := fastPolynomialReferenceCiphertext(params, input)
	require.NoError(t, fastEval.MulRelinElementQPrefixRows(input, input.El(), rows, t2))
	require.NoError(t, fastEval.AddQPrefixRows(t2, t2, t2, rows))
	require.NoError(t, fastEval.RescaleQPrefixRows(t2, rows, t2))
	rows, err = fastckks.QPrefixWidth(t2.Level())
	require.NoError(t, err)
	require.NoError(t, fastEval.AddScalarQPrefixRows(t2, -1, rows, t2))

	t3 := fastPolynomialReferenceCiphertext(params, input)
	require.NoError(t, fastEval.MulRelinElementQPrefixRows(t2, input.El(), rows, t3))
	require.NoError(t, fastEval.AddQPrefixRows(t3, t3, t3, rows))
	require.NoError(t, fastEval.RescaleQPrefixRows(t3, rows, t3))
	require.NoError(t, eval.workspace.subAligned(params, fastEval, t3, input))

	t5 := fastPolynomialReferenceCiphertext(params, input)
	require.NoError(t, fastEval.MulRelinElementQPrefixRows(t3, t2.El(), rows, t5))
	require.NoError(t, fastEval.AddQPrefixRows(t5, t5, t5, rows))
	require.NoError(t, fastEval.RescaleQPrefixRows(t5, rows, t5))
	require.NoError(t, eval.workspace.subAligned(params, fastEval, t5, input))

	require.Equal(t, t5.Level(), got.Level())
	require.True(t, t5.Scale.Equal(got.Scale), "reference scale=%v got=%v", t5.Scale.Float64(), got.Scale.Float64())
	require.Equal(t, t5.Degree(), got.Degree())
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < min(rows, got.Level()+1); limb++ {
			require.Equal(t, t5.Value[d].Coeffs[limb], got.Value[d].Coeffs[limb])
		}
	}
}

func fastPolynomialReferenceCiphertext(params ckks.Parameters, input *rlwe.Ciphertext) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, input.Level())
	*ct.MetaData = *input.MetaData
	ct.IsNTT = input.IsNTT
	ct.IsMontgomery = input.IsMontgomery
	return ct
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

func TestFastPolynomialPlanScaleOverrideKeepsInputAndFinalRescale(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 23)
	before := input.CopyNew()
	poly := fastPolynomialTestPoly(1)
	planScale := rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 29))
	got, err := eval.EvaluateWithPlanScale(input, poly, params.DefaultScale(), planScale)
	require.NoError(t, err)
	require.Equal(t, 1, got.Degree())
	require.Equal(t, input.Level()-1, got.Level())
	require.True(t, got.IsNTT)
	require.True(t, got.IsMontgomery)
	require.True(t, got.Scale.Cmp(params.DefaultScale()) < 0)
	for d := range input.Value {
		for limb := 0; limb <= input.Level(); limb++ {
			require.Equal(t, before.Value[d].Coeffs[limb], input.Value[d].Coeffs[limb])
		}
	}
}

func TestFastPolynomialFinalParentOneBitScalarGuardSelection(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 39, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56},
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	eval := NewFastEvaluator(params, nil)
	input := fastckks.NewCiphertext(params, 1, 12)
	input.IsNTT, input.IsMontgomery = true, true
	input.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	poly := bignum.NewPolynomial(bignum.Chebyshev, cosine.ApproximateCos(16, 30, float64(uint(1<<10)), 3), [2]float64{-16, 16})
	poly.IsOdd = false
	for i := range poly.Coeffs {
		if i&1 == 1 {
			poly.Coeffs[i] = nil
		}
	}
	planScale := rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 91))
	got, err := eval.EvaluateWithPlanScaleFinalParentOneBitScalarGuard(input, poly, rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60)), planScale)
	require.NoError(t, err)
	require.Equal(t, fastckks.MaintainedLimbCount(&params, input.Level()), eval.workspace.rows, "legacy guarded wrapper must preserve maintained-row authority")
	require.Equal(t, 1, got.Degree())
	require.True(t, got.IsNTT)
	require.True(t, got.IsMontgomery)
	require.True(t, eval.workspace.planScaleOverride)
	require.Equal(t, 1, eval.workspace.guardedOperations)
	require.Equal(t, eval.workspace.guardedPlanCount-1, eval.workspace.guardedPlanIndex)
	require.Equal(t, 2, eval.workspace.guardedScalarDegree)

	targetScale := rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	explicitEval := NewFastEvaluator(params, nil)
	explicit, err := explicitEval.EvaluateWithPlanScaleFinalParentOneBitScalarGuardQPrefixRows(input, poly, targetScale, planScale, fastckks.MaxQPrefixWidth)
	require.NoError(t, err)
	require.Equal(t, fastckks.MaxQPrefixWidth, explicitEval.workspace.rows)
	poisonedInput := input.CopyNew()
	q3 := params.Q()[3]
	for i, residue := range poisonedInput.Value[0].Coeffs[3] {
		poisonedInput.Value[0].Coeffs[3][i] = (residue + uint64(97+i)) % q3
	}
	explicitPoisonedEval := NewFastEvaluator(params, nil)
	explicitPoisoned, err := explicitPoisonedEval.EvaluateWithPlanScaleFinalParentOneBitScalarGuardQPrefixRows(poisonedInput, poly, targetScale, planScale, fastckks.MaxQPrefixWidth)
	require.NoError(t, err)
	require.NotEqual(t, explicit.Value[0].Coeffs[3], explicitPoisoned.Value[0].Coeffs[3], "explicit guarded wrapper must consume and produce q3")
}

func TestFastPolynomialIgnoresDormantResidues(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	poly := fastPolynomialTestPoly(3)
	inputA := fastPolynomialTestCiphertext(params, 19)
	inputB := inputA.CopyNew()
	inputRows, err := fastckks.QPrefixWidth(inputA.Level())
	require.NoError(t, err)
	for d := range inputB.Value {
		for limb := inputRows; limb <= inputB.Level(); limb++ {
			for i := range inputB.Value[d].Coeffs[limb] {
				inputB.Value[d].Coeffs[limb][i] ^= uint64(0x5a5a5a5a) + uint64(31*d+limb+i)
			}
		}
	}

	gotA, err := NewFastEvaluator(params, nil).Evaluate(inputA, poly, params.DefaultScale())
	require.NoError(t, err)
	gotB, err := NewFastEvaluator(params, nil).Evaluate(inputB, poly, params.DefaultScale())
	require.NoError(t, err)
	outputRows, err := fastckks.QPrefixWidth(gotA.Level())
	require.NoError(t, err)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < outputRows; limb++ {
			require.Equal(t, gotA.Value[d].Coeffs[limb], gotB.Value[d].Coeffs[limb])
		}
	}
}

func TestFastPolynomialQPrefixCopyAtLevel(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	src := fastPolynomialTestCiphertext(params, 29)
	for d := range src.Value {
		for limb := 2; limb <= src.Level(); limb++ {
			for i := range src.Value[d].Coeffs[limb] {
				src.Value[d].Coeffs[limb][i] = uint64(0xfeed0000 + 19*d + 7*limb + i)
			}
		}
	}
	srcBefore := src.CopyNew()
	dst := fastckks.NewCiphertext(params, 1, 1)
	require.NoError(t, copyQPrefixAtLevel(params, src, dst, 1, 2))

	require.Equal(t, 1, dst.Level())
	require.Equal(t, src.Scale, dst.Scale)
	require.Equal(t, src.IsNTT, dst.IsNTT)
	require.Equal(t, src.IsMontgomery, dst.IsMontgomery)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, src.Value[d].Coeffs[limb], dst.Value[d].Coeffs[limb])
		}
	}
	for d := range src.Value {
		for limb := range src.Value[d].Coeffs {
			require.Equal(t, srcBefore.Value[d].Coeffs[limb], src.Value[d].Coeffs[limb])
		}
	}
	require.Len(t, dst.Value[0].Coeffs, 2)
}

func TestFastPolynomialPublicResultCopiesQPrefixOnly(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 37)
	poly := fastPolynomialTestPoly(3)

	_, err := eval.Evaluate(input, poly, params.DefaultScale())
	require.NoError(t, err)
	workspaceResult := eval.workspace.babySteps[0].Value
	width, err := fastckks.QPrefixWidth(workspaceResult.Level())
	require.NoError(t, err)
	for d := range workspaceResult.Value {
		for limb := width; limb <= workspaceResult.Level(); limb++ {
			for i := range workspaceResult.Value[d].Coeffs[limb] {
				workspaceResult.Value[d].Coeffs[limb][i] = uint64(0xdead0000 + 101*d + 17*limb + i)
			}
		}
	}

	public := cloneQPrefixResult(params, workspaceResult, width)
	require.Equal(t, workspaceResult.Level(), public.Level())
	require.Equal(t, workspaceResult.Scale, public.Scale)
	require.Equal(t, workspaceResult.IsNTT, public.IsNTT)
	require.Equal(t, workspaceResult.IsMontgomery, public.IsMontgomery)
	for d := range workspaceResult.Value {
		for limb := 0; limb < width; limb++ {
			require.Equal(t, workspaceResult.Value[d].Coeffs[limb], public.Value[d].Coeffs[limb])
		}
		for limb := width; limb <= public.Level(); limb++ {
			require.Empty(t, public.Value[d].Coeffs[limb])
		}
	}
}

func TestFastPolynomialQ3PoisonPropagatesAcrossBabyGiantAndFinalPS(t *testing.T) {
	params := fastPolynomialQ012TestParameters(t)
	require.Equal(t, 3, fastckks.MaintainedLimbCount(&params, params.MaxLevel()), "fixture must exercise legacy q012 authority; Q0..Q2=%v", params.Q()[:3])
	const rows = fastckks.MaxQPrefixWidth
	poisonQ3 := func(ct *rlwe.Ciphertext) *rlwe.Ciphertext {
		poisoned := ct.CopyNew()
		q3 := params.Q()[3]
		for i, residue := range poisoned.Value[0].Coeffs[3] {
			poisoned.Value[0].Coeffs[3][i] = (residue + uint64(97+i)) % q3
		}
		return poisoned
	}
	assertQ012SameQ3Different := func(clean, poisoned *rlwe.Ciphertext, context string) {
		t.Helper()
		require.Equal(t, clean.Level(), poisoned.Level(), context)
		q3Changed := false
		for component := range clean.Value {
			for row := 0; row < 3; row++ {
				require.Equal(t, clean.Value[component].Coeffs[row], poisoned.Value[component].Coeffs[row], "%s q%d component=%d", context, row, component)
			}
			for i, residue := range clean.Value[component].Coeffs[3] {
				if residue != poisoned.Value[component].Coeffs[3][i] {
					q3Changed = true
					break
				}
			}
		}
		require.True(t, q3Changed, "%s must consume and freshly produce q3", context)
	}

	// Baby-step accumulation has no Rescale, so q3-only input poison must
	// remain rowwise: q012 are identical and q3 changes.
	baby1 := fastPolynomialTestCiphertext(params, 101)
	poisonedBaby1 := poisonQ3(baby1)
	baby2 := fastPolynomialTestCiphertext(params, 151)
	poly := commonpolynomial.NewPolynomial(fastPolynomialTestPoly(2))
	poly.Level, poly.Scale = baby1.Level(), baby1.Scale
	cleanBabyEval, poisonedBabyEval := NewFastEvaluator(params, nil), NewFastEvaluator(params, nil)
	cleanBabyEval.workspace.reset(params, baby1, rows)
	poisonedBabyEval.workspace.reset(params, poisonedBaby1, rows)
	cleanBaby, err := cleanBabyEval.workspace.evaluateBabyStep(params, cleanBabyEval.Evaluator, poly, map[int]*rlwe.Ciphertext{1: baby1, 2: baby2}, 0, false)
	require.NoError(t, err)
	poisonedBaby, err := poisonedBabyEval.workspace.evaluateBabyStep(params, poisonedBabyEval.Evaluator, poly, map[int]*rlwe.Ciphertext{1: poisonedBaby1, 2: baby2}, 0, false)
	require.NoError(t, err)
	assertQ012SameQ3Different(cleanBaby.Value, poisonedBaby.Value, "baby-step")

	// In a giant step the q3-only poison is placed in xpow, after the
	// independently identical b Rescale. The following multiply/add must
	// update q3 without cross-row contamination of q012.
	newAtLevel := func(offset uint64, level int) *rlwe.Ciphertext {
		ct := fastPolynomialTestCiphertext(params, offset)
		fastckks.Resize(ct, ct.Degree(), level, params.N())
		return ct
	}
	makeGiantOperands := func(poisoned bool) (a, b, xpow *rlwe.Ciphertext) {
		b = newAtLevel(181, 4)
		xpow = newAtLevel(211, 3)
		a = newAtLevel(241, 3)
		b.Scale = rlwe.NewScale(1 << 40)
		xpow.Scale = rlwe.NewScale(1 << 30)
		a.Scale = b.Scale.Div(rlwe.NewScale(params.Q()[4])).Mul(xpow.Scale)
		if poisoned {
			xpow = poisonQ3(xpow)
		}
		return
	}
	cleanA, cleanB, cleanXpow := makeGiantOperands(false)
	poisonedA, poisonedB, poisonedXpow := makeGiantOperands(true)
	cleanGiantEval, poisonedGiantEval := NewFastEvaluator(params, nil), NewFastEvaluator(params, nil)
	cleanGiantEval.workspace.reset(params, cleanB, rows)
	poisonedGiantEval.workspace.reset(params, poisonedB, rows)
	require.NoError(t, cleanGiantEval.workspace.evaluateMonomial(params, cleanGiantEval.Evaluator, cleanA, cleanB, cleanXpow))
	require.NoError(t, poisonedGiantEval.workspace.evaluateMonomial(params, poisonedGiantEval.Evaluator, poisonedA, poisonedB, poisonedXpow))
	assertQ012SameQ3Different(cleanB, poisonedB, "giant-step")

	// Public PS output must not inherit a stale q3 row when the input q3
	// history differs.
	finalInput := fastPolynomialTestCiphertext(params, 271)
	polynomial := fastPolynomialTestPoly(3)
	planScale := rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 29))
	for _, testCase := range []struct {
		name     string
		legacy   func(*FastEvaluator, *rlwe.Ciphertext) (*rlwe.Ciphertext, error)
		explicit func(*FastEvaluator, *rlwe.Ciphertext) (*rlwe.Ciphertext, error)
	}{
		{
			name: "evaluate",
			legacy: func(eval *FastEvaluator, input *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
				return eval.Evaluate(input, polynomial, params.DefaultScale())
			},
			explicit: func(eval *FastEvaluator, input *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
				return eval.EvaluateQPrefixRows(input, polynomial, params.DefaultScale(), rows)
			},
		},
		{
			name: "plan-scale",
			legacy: func(eval *FastEvaluator, input *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
				return eval.EvaluateWithPlanScale(input, fastPolynomialTestPoly(1), params.DefaultScale(), planScale)
			},
			explicit: func(eval *FastEvaluator, input *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
				return eval.EvaluateWithPlanScaleQPrefixRows(input, fastPolynomialTestPoly(1), params.DefaultScale(), planScale, rows)
			},
		},
	} {
		t.Run(testCase.name, func(t *testing.T) {
			legacyFinal, err := testCase.legacy(NewFastEvaluator(params, nil), finalInput)
			require.NoError(t, err)
			legacyPoisoned, err := testCase.legacy(NewFastEvaluator(params, nil), poisonQ3(finalInput))
			require.NoError(t, err)
			require.Equal(t, legacyFinal.Level(), legacyPoisoned.Level())
			for d := 0; d <= 1; d++ {
				for row := 0; row < 3; row++ {
					require.Equal(t, legacyFinal.Value[d].Coeffs[row], legacyPoisoned.Value[d].Coeffs[row], "legacy wrapper q%d component=%d", row, d)
				}
			}

			explicitFinal, err := testCase.explicit(NewFastEvaluator(params, nil), finalInput)
			require.NoError(t, err)
			explicitPoisoned, err := testCase.explicit(NewFastEvaluator(params, nil), poisonQ3(finalInput))
			require.NoError(t, err)
			require.Equal(t, explicitFinal.Level(), explicitPoisoned.Level())
			require.NotEqual(t, explicitFinal.Value[0].Coeffs[3], explicitPoisoned.Value[0].Coeffs[3], "explicit-row wrapper must consume and produce q3")
		})
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

func TestFastPolynomialPatersonStockmeyerNumericalOracle(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            6,
		LogQ:            []int{55, 39, 39, 39, 39, 39, 39, 39},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	encoder := ckks.NewEncoder(params)
	inputLevel := params.MaxLevel() - 1
	pt := ckks.NewPlaintext(params, inputLevel)
	pt.IsNTT = true
	pt.IsMontgomery = true
	values := make([]complex128, pt.Slots())
	for i := range values {
		values[i] = complex(-0.02+0.001*float64(i), 0.004-0.0003*float64(i))
	}
	require.NoError(t, encoder.Encode(values, pt))

	input := ckks.NewCiphertext(params, 1, inputLevel)
	*input.MetaData = *pt.MetaData
	input.Value[0].Copy(pt.Value)
	input.Value[1].Zero()

	coeffs := make([]*bignum.Complex, 31)
	for i := range coeffs {
		coeffs[i] = bignum.ToComplex(0.08*math.Cos(float64(i+1))/(float64(i)+1), params.EncodingPrecision())
	}
	poly := bignum.NewPolynomial(bignum.Chebyshev, coeffs, [2]float64{-1, 1})

	fastEval := NewFastEvaluator(params, nil)
	got, err := fastEval.Evaluate(input, poly, params.DefaultScale())
	require.NoError(t, err)
	require.Equal(t, 1, got.Level())
	require.True(t, got.Scale.Equal(params.DefaultScale()))

	decoded := decodeFastPolynomialOutput(t, params, encoder, got)
	want := make([]complex128, len(values))
	for i := range values {
		x := bignum.ToComplex(values[i], params.EncodingPrecision())
		wantComplex := poly.Clone()
		want[i] = complexFloat64(wantComplex.Evaluate(x))
	}
	for i := range want {
		// The q1=39-bit test chain and five rescaling stages leave a
		// deliberately modest CKKS tolerance for this execution oracle.
		require.InDelta(t, real(want[i]), real(decoded[i]), 1e-2)
		require.InDelta(t, imag(want[i]), imag(decoded[i]), 1e-2)
	}
}

func TestFastPolynomialFormalScaleT2CapacitySafe(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 39, 39},
		LogDefaultScale: 60,
	})
	require.NoError(t, err)
	encoder := ckks.NewEncoder(params)
	inputLevel := params.MaxLevel()
	pt := ckks.NewPlaintext(params, inputLevel)
	pt.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	pt.IsNTT = true
	pt.IsMontgomery = true
	values := make([]complex128, pt.Slots())
	for i := range values {
		values[i] = complex(float64((i%5)-2)*1e-5, float64((i%3)-1)*1e-5)
	}
	require.NoError(t, encoder.Encode(values, pt))

	input := ckks.NewCiphertext(params, 1, inputLevel)
	*input.MetaData = *pt.MetaData
	input.Value[0].Copy(pt.Value)
	input.Value[1].Zero()
	input.IsNTT = pt.IsNTT
	input.IsMontgomery = pt.IsMontgomery

	poly := bignum.NewPolynomial(bignum.Chebyshev, []*bignum.Complex{
		bignum.ToComplex(0, params.EncodingPrecision()),
		bignum.ToComplex(0, params.EncodingPrecision()),
		bignum.ToComplex(1, params.EncodingPrecision()),
	}, [2]float64{-1, 1})
	eval := NewFastEvaluator(params, nil)
	commonPoly := commonpolynomial.NewPolynomial(poly)
	rows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)
	eval.workspace.reset(params, input, rows)
	require.NoError(t, eval.workspace.generatePowers(params, eval.Evaluator, poly, commonPoly))
	got := eval.workspace.powers[2]
	require.Equal(t, inputLevel-1, got.Level())

	decoded := decodeFastPolynomialOutput(t, params, encoder, got)
	for i := range values {
		want := 2*values[i]*values[i] - 1
		require.InDelta(t, real(want), real(decoded[i]), 1e-2)
		require.InDelta(t, imag(want), imag(decoded[i]), 1e-2)
	}
}

func TestFastPolynomialFormalScaleT3CapacitySafe(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 39, 39, 39},
		LogDefaultScale: 60,
	})
	require.NoError(t, err)
	encoder := ckks.NewEncoder(params)
	inputLevel := params.MaxLevel()
	pt := ckks.NewPlaintext(params, inputLevel)
	pt.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	pt.IsNTT = true
	pt.IsMontgomery = true
	values := make([]complex128, pt.Slots())
	for i := range values {
		values[i] = complex(float64((i%5)-2)/8192, float64((i%3)-1)/8192)
	}
	require.NoError(t, encoder.Encode(values, pt))

	input := ckks.NewCiphertext(params, 1, inputLevel)
	*input.MetaData = *pt.MetaData
	input.Value[0].Copy(pt.Value)
	input.Value[1].Zero()
	input.IsNTT = pt.IsNTT
	input.IsMontgomery = pt.IsMontgomery

	eval := NewFastEvaluator(params, nil)
	rows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)
	eval.workspace.reset(params, input, rows)
	pb := fastPowerBasis{basis: bignum.Chebyshev, values: eval.workspace.powers, workspace: &eval.workspace, params: params, eval: eval.Evaluator}
	require.NoError(t, pb.genPower(3, false))
	got := eval.workspace.powers[3]
	require.Equal(t, inputLevel-2, got.Level())

	decoded := decodeFastPolynomialOutput(t, params, encoder, got)
	for i := range values {
		want := 4*values[i]*values[i]*values[i] - 3*values[i]
		require.InDelta(t, real(want), real(decoded[i]), 1e-2)
		require.InDelta(t, imag(want), imag(decoded[i]), 1e-2)
	}
}

func TestFastPolynomialBalancedLazyMetadata(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            6,
		LogQ:            []int{55, 39, 39, 39, 39, 39},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 47)
	input.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	rows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)
	eval.workspace.reset(params, input, rows)
	pb := fastPowerBasis{basis: bignum.Chebyshev, values: eval.workspace.powers, workspace: &eval.workspace, params: params, eval: eval.Evaluator}
	require.NoError(t, pb.genPower(3, true))
	got := eval.workspace.powers[3]
	require.Equal(t, 2, got.Degree())
	require.Equal(t, input.Level()-2, got.Level())
	wantScale := input.Scale.Mul(input.Scale).Mul(input.Scale)
	wantScale = wantScale.Div(rlwe.NewScale(params.Q()[input.Level()]))
	wantScale = wantScale.Div(rlwe.NewScale(params.Q()[input.Level()-1]))
	require.True(t, got.Scale.Equal(wantScale), "unexpected lazy T3 scale: got=%v want=%v", got.Scale.Float64(), wantScale.Float64())
}

func TestBalancedFactorPairDeterministic(t *testing.T) {
	formalQ := uint64(1152921504607191041)
	factors, err := balancedFactorPair(formalQ)
	require.NoError(t, err)
	require.Equal(t, uint64(1073741824), factors.left)
	require.Equal(t, uint64(1073741824), factors.right)

	for _, q := range []uint64{(uint64(1) << 60) - 33, (uint64(1) << 60) + 344065, (uint64(1) << 61) - 12345} {
		factors, err := balancedFactorPair(q)
		require.NoError(t, err)
		product := new(big.Int).Mul(new(big.Int).SetUint64(factors.left), new(big.Int).SetUint64(factors.right))
		require.True(t, rlwe.NewScale(product).InDelta(rlwe.NewScale(q), balancedScaleToleranceBits))
	}
}

func TestBalancedScheduleBranchSelection(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	left := fastPolynomialTestCiphertext(params, 7)
	right := fastPolynomialTestCiphertext(params, 11)
	left.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	right.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	pb := fastPowerBasis{params: params, eval: eval.Evaluator}
	high, err := pb.balancedScheduleFor(left, right, left.Level())
	require.NoError(t, err)
	require.True(t, high.balanced)
	require.GreaterOrEqual(t, high.leftScale.Cmp(rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), balancedMinScaleBits))), 0)

	left.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 30))
	right.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 30))
	low, err := pb.balancedScheduleFor(left, right, left.Level())
	require.NoError(t, err)
	require.False(t, low.balanced)
}

func TestPostProductQ012ScheduleSelection(t *testing.T) {
	q012Params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		Q:               []uint64{72057594037616641, 549755731969, 549756026881},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	q012 := fastPowerBasis{basis: bignum.Chebyshev, params: q012Params}
	require.False(t, q012.postProductLogN13Schedule(1), "the validated schedule must not be selected before its level floor")
	require.True(t, q012.postProductLogN13Schedule(2))

	legacyParams, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 39, 40},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	legacy := fastPowerBasis{basis: bignum.Chebyshev, params: legacyParams}
	require.False(t, legacy.postProductLogN13Schedule(2), "non-P93 profiles must retain the fallback schedule")

	monomial := fastPowerBasis{basis: bignum.Monomial, params: q012Params}
	require.False(t, monomial.postProductLogN13Schedule(2), "the post-product recurrence is Chebyshev-specific")
}

func TestPostProductQ012PowerContract(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{56, 39, 40, 40},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 17)
	rows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)
	eval.workspace.reset(params, input, rows)
	pb := fastPowerBasis{basis: bignum.Chebyshev, values: eval.workspace.powers, workspace: &eval.workspace, params: params, eval: eval.Evaluator}
	require.NoError(t, pb.genPower(2, false))
	got := eval.workspace.powers[2]

	reference := fastPolynomialReferenceCiphertext(params, input)
	require.NoError(t, eval.Evaluator.MulRelinElementQPrefixRows(input, input.El(), rows, reference))
	require.NoError(t, eval.Evaluator.AddQPrefixRows(reference, reference, reference, rows))
	require.NoError(t, eval.Evaluator.AddScalarQPrefixRows(reference, -1, rows, reference))
	require.NoError(t, eval.Evaluator.RescaleQPrefixRows(reference, rows, reference))

	require.Equal(t, input.Level()-1, got.Level(), "the repaired recurrence consumes exactly one level")
	require.Equal(t, 1, got.Degree())
	require.True(t, got.Scale.Equal(reference.Scale), "post-product scale=%v reference=%v", got.Scale.Float64(), reference.Scale.Float64())
	rows, err = fastckks.QPrefixWidth(got.Level())
	require.NoError(t, err)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < rows; limb++ {
			require.Equal(t, reference.Value[d].Coeffs[limb], got.Value[d].Coeffs[limb], "component=%d limb=%d", d, limb)
		}
	}
}

func TestBalancedPowerSourceImmutabilityAndScratchOwnership(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            6,
		LogQ:            []int{55, 39, 39, 39, 39, 39},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 7)
	input.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	rows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)
	eval.workspace.reset(params, input, rows)
	pb := fastPowerBasis{basis: bignum.Chebyshev, values: eval.workspace.powers, workspace: &eval.workspace, params: params, eval: eval.Evaluator}
	require.NoError(t, pb.genPower(2, false))
	leftBefore := eval.workspace.powers[1].CopyNew()
	rightBefore := eval.workspace.powers[2].CopyNew()
	require.NoError(t, pb.genPower(3, false))
	for d := range leftBefore.Value {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, leftBefore.Value[d].Coeffs[limb], eval.workspace.powers[1].Value[d].Coeffs[limb])
			require.Equal(t, rightBefore.Value[d].Coeffs[limb], eval.workspace.powers[2].Value[d].Coeffs[limb])
		}
	}
	require.NotSame(t, eval.workspace.balancedLeft, eval.workspace.balancedRight)
	require.NotSame(t, eval.workspace.balancedLeft, eval.workspace.powers[1])
	require.NotSame(t, eval.workspace.balancedRight, eval.workspace.powers[2])
	require.Equal(t, 1, eval.workspace.powers[3].Degree())
	require.Equal(t, input.Level()-2, eval.workspace.powers[3].Level())
	require.True(t, eval.workspace.powers[3].Scale.Equal(input.Scale.Mul(input.Scale).Mul(input.Scale).Div(rlwe.NewScale(params.Q()[input.Level()])).Div(rlwe.NewScale(params.Q()[input.Level()-1]))))
}

func TestBalancedPowerIgnoresDormantResidues(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            6,
		LogQ:            []int{55, 39, 39, 39, 39, 39},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	inputA := fastPolynomialTestCiphertext(params, 17)
	inputA.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	inputB := inputA.CopyNew()
	inputRows, err := fastckks.QPrefixWidth(inputA.Level())
	require.NoError(t, err)
	for d := range inputB.Value {
		for limb := inputRows; limb <= inputB.Level(); limb++ {
			for i := range inputB.Value[d].Coeffs[limb] {
				inputB.Value[d].Coeffs[limb][i] ^= uint64(0x9e3779b9) + uint64(13*d+limb+i)
			}
		}
	}
	run := func(input *rlwe.Ciphertext) *rlwe.Ciphertext {
		eval := NewFastEvaluator(params, nil)
		rows, err := fastckks.QPrefixWidth(input.Level())
		require.NoError(t, err)
		eval.workspace.reset(params, input, rows)
		pb := fastPowerBasis{basis: bignum.Chebyshev, values: eval.workspace.powers, workspace: &eval.workspace, params: params, eval: eval.Evaluator}
		require.NoError(t, pb.genPower(3, false))
		return eval.workspace.powers[3]
	}
	gotA := run(inputA)
	gotB := run(inputB)
	outputRows, err := fastckks.QPrefixWidth(gotA.Level())
	require.NoError(t, err)
	for d := range gotA.Value {
		for limb := 0; limb < outputRows; limb++ {
			require.Equal(t, gotA.Value[d].Coeffs[limb], gotB.Value[d].Coeffs[limb])
		}
	}
}

func decodeFastPolynomialOutput(t *testing.T, params ckks.Parameters, encoder *ckks.Encoder, ct *rlwe.Ciphertext) []complex128 {
	t.Helper()
	pt := ckks.NewPlaintext(params, ct.Level())
	*pt.MetaData = *ct.MetaData
	pt.Value.Copy(ct.Value[0])
	if pt.IsMontgomery {
		params.RingQ().AtLevel(pt.Level()).INTT(pt.Value, pt.Value)
		params.RingQ().AtLevel(pt.Level()).IMForm(pt.Value, pt.Value)
		pt.IsNTT = false
		pt.IsMontgomery = false
	}
	values := make([]complex128, pt.Slots())
	require.NoError(t, encoder.Decode(pt, values))
	return values
}

func complexFloat64(value *bignum.Complex) complex128 {
	realValue, _ := value[0].Float64()
	imagValue, _ := value[1].Float64()
	return complex(realValue, imagValue)
}

func TestFastPolynomialEvaluatorValidation(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 1)
	poly := fastPolynomialTestPoly(3)
	_, err := eval.Evaluate(input, bignum.Polynomial{}, params.DefaultScale())
	require.Error(t, err)
	_, err = eval.EvaluateQPrefixRows(input, poly, params.DefaultScale(), 0)
	require.Error(t, err, "an invalid explicit row count must not fall back to legacy authority")
	_, err = eval.EvaluateQPrefixRows(input, poly, params.DefaultScale(), -1)
	require.Error(t, err, "negative explicit rows must not select legacy authority")

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
