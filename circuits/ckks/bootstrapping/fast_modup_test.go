package bootstrapping

import (
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func fastModUpParameters(t testing.TB, logN int) (Parameters, ckks.Parameters) {
	t.Helper()
	logQ := []int{55, 39, 50, 50}
	if logN > 4 {
		logQ = []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56}
		if logN == 16 {
			logQ[0], logQ[1] = 54, 38
		}
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            logQ,
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	return Parameters{
		BootstrappingParameters: params,
		Mod1ParametersLiteral: mod1.ParametersLiteral{
			LogScale:        30,
			LogMessageRatio: 2,
		},
	}, params
}

func newFastModUpInput(params ckks.Parameters, values []int64, scale rlwe.Scale, withBacking bool) *rlwe.Ciphertext {
	level := 0
	if withBacking {
		level = params.MaxLevel()
	}
	ct := ckks.NewCiphertext(params, 1, level)
	ct.IsNTT = true
	ct.Scale = scale
	q0 := params.RingQ().SubRings[0].Modulus
	for component := range ct.Value {
		for j := range ct.Value[component].Coeffs[0] {
			value := values[(j+component)%len(values)]
			ct.Value[component].Coeffs[0][j] = new(big.Int).Mod(big.NewInt(value), new(big.Int).SetUint64(q0)).Uint64()
		}
		params.RingQ().SubRings[0].NTT(ct.Value[component].Coeffs[0], ct.Value[component].Coeffs[0])
		if withBacking {
			for limb := 1; limb <= level; limb++ {
				for j := range ct.Value[component].Coeffs[limb] {
					ct.Value[component].Coeffs[limb][j] = uint64(0x9e3779b9) + uint64(limb*params.N()+j+component)
				}
			}
		}
	}
	if withBacking {
		ct.Resize(ct.Degree(), 0)
	}
	return ct
}

func standardModUpBasisReference(params Parameters, ct *rlwe.Ciphertext) *rlwe.Ciphertext {
	ringQ := params.BootstrappingParameters.RingQ()
	for component := range ct.Value {
		ringQ.AtLevel(0).INTT(ct.Value[component], ct.Value[component])
	}
	ct.Resize(ct.Degree(), params.BootstrappingParameters.MaxLevel())
	Q := ringQ.ModuliChain()
	q0 := Q[0]
	BRCQ := ringQ.BRedConstants()
	for component := range ct.Value {
		for j := range ct.Value[component].Coeffs[0] {
			coeff := ct.Value[component].Coeffs[0][j]
			pos, neg := uint64(1), uint64(0)
			if coeff >= q0>>1 {
				coeff = q0 - coeff
				pos, neg = 0, 1
			}
			for limb := 1; limb <= params.BootstrappingParameters.MaxLevel(); limb++ {
				tmp := ring.BRedAdd(coeff, Q[limb], BRCQ[limb])
				ct.Value[component].Coeffs[limb][j] = tmp*pos + (Q[limb]-tmp)*neg
			}
		}
		ringQ.NTT(ct.Value[component], ct.Value[component])
	}
	mod1Parameters := mod1.Parameters{
		LogDefaultScale: params.Mod1ParametersLiteral.LogScale,
		LogMessageRatio: params.Mod1ParametersLiteral.LogMessageRatio,
	}
	if scale := (mod1Parameters.ScalingFactor().Float64() / mod1Parameters.MessageRatio()) / ct.Scale.Float64(); scale > 1 {
		scalar := uint64(math.Round(scale))
		ringQ.MulScalar(ct.Value[0], scalar, ct.Value[0])
		ringQ.MulScalar(ct.Value[1], scalar, ct.Value[1])
		ct.Scale = ct.Scale.Mul(rlwe.NewScale(scale))
	}
	return ct
}

func requireFastModUpBasisMatches(t *testing.T, want, got *rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, want.Level(), got.Level())
	require.True(t, want.Scale.Equal(got.Scale))
	require.Equal(t, want.IsNTT, got.IsNTT)
	require.Equal(t, want.IsMontgomery, got.IsMontgomery)
	for component := range got.Value {
		require.Equal(t, want.Value[component].Coeffs[0], got.Value[component].Coeffs[0], "component %d q0", component)
		require.Equal(t, want.Value[component].Coeffs[1], got.Value[component].Coeffs[1], "component %d q1", component)
	}
}

func TestFastModUpBasisMatchesStandardCenteredLiftAndScale(t *testing.T) {
	params, ckksParams := fastModUpParameters(t, 4)
	q0 := ckksParams.RingQ().SubRings[0].Modulus
	values := []int64{
		int64(q0>>1) - 1,
		int64(q0 >> 1),
		-int64(q0>>1) + 1,
		-int64(q0 >> 1),
		12345,
		-67890,
	}
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	input := newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<20), true)
	want := standardModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<20), false))
	got, err := fastEval.modUpBasis(input)
	require.NoError(t, err)
	requireFastModUpBasisMatches(t, want, got)
	require.Equal(t, params.BootstrappingParameters.MaxLevel(), got.Level())
}

func TestFastModUpBasisNoScaleUpAndPoisonedDormantRows(t *testing.T) {
	params, ckksParams := fastModUpParameters(t, 4)
	values := []int64{1, -2, 12345, -67890}
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	input := ckks.NewCiphertext(ckksParams, 1, ckksParams.MaxLevel())
	input.IsNTT = true
	input.Scale = rlwe.NewScale(1 << 40)
	q0 := ckksParams.RingQ().SubRings[0].Modulus
	rowPointers := make([][]*uint64, 2)
	for component := range input.Value {
		for j := range input.Value[component].Coeffs[0] {
			value := values[(j+component)%len(values)]
			input.Value[component].Coeffs[0][j] = new(big.Int).Mod(big.NewInt(value), new(big.Int).SetUint64(q0)).Uint64()
		}
		ckksParams.RingQ().SubRings[0].NTT(input.Value[component].Coeffs[0], input.Value[component].Coeffs[0])
		rowPointers[component] = make([]*uint64, ckksParams.MaxLevel()-1)
		for limb := 2; limb <= ckksParams.MaxLevel(); limb++ {
			for j := range input.Value[component].Coeffs[limb] {
				input.Value[component].Coeffs[limb][j] = uint64(0x9e3779b9) + uint64(limb*ckksParams.N()+j+component)
			}
			rowPointers[component][limb-2] = &input.Value[component].Coeffs[limb][0]
		}
	}
	input.Resize(input.Degree(), 0)
	want := standardModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<40), false))
	got, err := fastEval.modUpBasis(input)
	require.NoError(t, err)
	requireFastModUpBasisMatches(t, want, got)
	require.Equal(t, rlwe.NewScale(1<<40), got.Scale)
	for component := range got.Value {
		for limb := 2; limb <= ckksParams.MaxLevel(); limb++ {
			require.Equal(t, rowPointers[component][limb-2], &got.Value[component].Coeffs[limb][0])
			poison := uint64(0x9e3779b9) + uint64(limb*ckksParams.N()+component)
			require.Equal(t, poison, got.Value[component].Coeffs[limb][0])
		}
	}
}

func TestFastModUpBasisFreshLevelZeroStorage(t *testing.T) {
	params, ckksParams := fastModUpParameters(t, 4)
	values := []int64{7, -11, 1234, -5678}
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	got, err := fastEval.modUpBasis(newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<40), false))
	require.NoError(t, err)
	want := standardModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<40), false))
	requireFastModUpBasisMatches(t, want, got)
}

func TestFastModUpBasisValidation(t *testing.T) {
	params, ckksParams := fastModUpParameters(t, 4)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	_, err = fastEval.modUpBasis(nil)
	require.Error(t, err)

	wrongDegree := ckks.NewCiphertext(ckksParams, 2, 0)
	wrongDegree.IsNTT = true
	_, err = fastEval.modUpBasis(wrongDegree)
	require.Error(t, err)

	wrongLevel := ckks.NewCiphertext(ckksParams, 1, 1)
	wrongLevel.IsNTT = true
	_, err = fastEval.modUpBasis(wrongLevel)
	require.Error(t, err)

	nonNTT := newFastModUpInput(ckksParams, []int64{1, 2}, rlwe.NewScale(1), false)
	nonNTT.IsNTT = false
	_, err = fastEval.modUpBasis(nonNTT)
	require.Error(t, err)

	montgomery := newFastModUpInput(ckksParams, []int64{1, 2}, rlwe.NewScale(1), false)
	montgomery.IsMontgomery = true
	_, err = fastEval.modUpBasis(montgomery)
	require.Error(t, err)

	ciParams, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     4,
		LogQ:     []int{55, 39, 50, 50},
		RingType: ring.ConjugateInvariant,
	})
	require.NoError(t, err)
	_, err = NewFastEvaluator(Parameters{BootstrappingParameters: ciParams})
	require.Error(t, err)
}

func resetFastModUpInput(ct, template *rlwe.Ciphertext, scale rlwe.Scale) {
	for component := range ct.Value {
		ct.Value[component].Coeffs = ct.Value[component].Coeffs[:1]
		copy(ct.Value[component].Coeffs[0], template.Value[component].Coeffs[0])
	}
	ct.Scale = scale
}

func BenchmarkFastModUpBasisLogN13(b *testing.B) {
	benchmarkFastModUpBasis(b, 13, true)
}

func BenchmarkStandardModUpBasisLogN13(b *testing.B) {
	benchmarkFastModUpBasis(b, 13, false)
}

func BenchmarkFastModUpBasisLogN16(b *testing.B) {
	benchmarkFastModUpBasis(b, 16, true)
}

func BenchmarkStandardModUpBasisLogN16(b *testing.B) {
	benchmarkFastModUpBasis(b, 16, false)
}

func benchmarkFastModUpBasis(b *testing.B, logN int, fastPath bool) {
	params, ckksParams := fastModUpParameters(b, logN)
	values := []int64{1 << 40, -(1 << 40) + 17, 1 << 39, -123456789}
	scale := rlwe.NewScale(1 << 40)
	template := newFastModUpInput(ckksParams, values, scale, true)
	input := newFastModUpInput(ckksParams, values, scale, true)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(b, err)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		resetFastModUpInput(input, template, scale)
		b.StartTimer()
		if fastPath {
			if _, err := fastEval.modUpBasis(input); err != nil {
				b.Fatal(err)
			}
		} else {
			standardModUpBasisReference(params, input)
		}
	}
}

func BenchmarkFastModUpBasisColdLogN13(b *testing.B) {
	params, ckksParams := fastModUpParameters(b, 13)
	values := []int64{1 << 40, -(1 << 40) + 17, 1 << 39, -123456789}
	scale := rlwe.NewScale(1 << 40)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(b, err)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		input := newFastModUpInput(ckksParams, values, scale, false)
		b.StartTimer()
		if _, err := fastEval.modUpBasis(input); err != nil {
			b.Fatal(err)
		}
	}
}
