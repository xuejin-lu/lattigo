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
	require.NoError(t, fastEval.Rescale(reference, reference))
	require.NoError(t, fastEval.Add(reference, -1, reference))

	require.Equal(t, reference.Level(), got.Level())
	require.True(t, reference.Scale.Equal(got.Scale), "reference scale=%v got=%v", reference.Scale.Float64(), got.Scale.Float64())
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
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

	_, err = eval.Evaluate(input, poly, params.DefaultScale())
	require.NoError(t, err)
	got := eval.workspace.powers[5]
	fastEval := eval.Evaluator

	t2 := fastPolynomialReferenceCiphertext(params, input)
	require.NoError(t, fastEval.MulRelin(input, input, t2))
	require.NoError(t, fastEval.Add(t2, t2, t2))
	require.NoError(t, fastEval.Rescale(t2, t2))
	require.NoError(t, fastEval.Add(t2, -1, t2))

	t3 := fastPolynomialReferenceCiphertext(params, input)
	require.NoError(t, fastEval.MulRelin(t2, input, t3))
	require.NoError(t, fastEval.Add(t3, t3, t3))
	require.NoError(t, fastEval.Rescale(t3, t3))
	require.NoError(t, eval.workspace.subAligned(params, fastEval, t3, input))

	t5 := fastPolynomialReferenceCiphertext(params, input)
	require.NoError(t, fastEval.MulRelin(t3, t2, t5))
	require.NoError(t, fastEval.Add(t5, t5, t5))
	require.NoError(t, fastEval.Rescale(t5, t5))
	require.NoError(t, eval.workspace.subAligned(params, fastEval, t5, input))

	require.Equal(t, t5.Level(), got.Level())
	require.True(t, t5.Scale.Equal(got.Scale), "reference scale=%v got=%v", t5.Scale.Float64(), got.Scale.Float64())
	require.Equal(t, t5.Degree(), got.Degree())
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < 2; limb++ {
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

func TestFastPolynomialMaintainedCopyAtLevel(t *testing.T) {
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
	require.NoError(t, copyMaintainedAtLevel(params, src, dst, 1))

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

func TestFastPolynomialPublicResultCopiesMaintainedResiduesOnly(t *testing.T) {
	params := fastPolynomialTestParameters(t)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 37)
	poly := fastPolynomialTestPoly(3)

	_, err := eval.Evaluate(input, poly, params.DefaultScale())
	require.NoError(t, err)
	workspaceResult := eval.workspace.babySteps[0].Value
	for d := range workspaceResult.Value {
		for limb := 2; limb <= workspaceResult.Level(); limb++ {
			for i := range workspaceResult.Value[d].Coeffs[limb] {
				workspaceResult.Value[d].Coeffs[limb][i] = uint64(0xdead0000 + 101*d + 17*limb + i)
			}
		}
	}

	public := cloneMaintainedResult(params, workspaceResult)
	require.Equal(t, workspaceResult.Level(), public.Level())
	require.Equal(t, workspaceResult.Scale, public.Scale)
	require.Equal(t, workspaceResult.IsNTT, public.IsNTT)
	require.Equal(t, workspaceResult.IsMontgomery, public.IsMontgomery)
	for d := range workspaceResult.Value {
		for limb := 0; limb < 2; limb++ {
			require.Equal(t, workspaceResult.Value[d].Coeffs[limb], public.Value[d].Coeffs[limb])
		}
		for limb := 2; limb <= public.Level(); limb++ {
			require.Empty(t, public.Value[d].Coeffs[limb])
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
	eval.workspace.reset(params, input)
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
	eval.workspace.reset(params, input)
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

func TestFastPolynomialPreRescaleLazyMetadata(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            6,
		LogQ:            []int{55, 39, 39, 39, 39, 39},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	eval := NewFastEvaluator(params, nil)
	input := fastPolynomialTestCiphertext(params, 47)
	input.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 60))
	eval.workspace.reset(params, input)
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
