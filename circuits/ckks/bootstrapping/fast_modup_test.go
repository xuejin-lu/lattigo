package bootstrapping

import (
	"fmt"
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
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

func fastModUpQ012Parameters(t testing.TB) (Parameters, ckks.Parameters) {
	t.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(45, 2*2*16)
	q3, err := generator.NextAlternatingPrimes(1)
	require.NoError(t, err)
	ckksParams, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		Q:               []uint64{72057594037616641, 549755731969, 549756026881, q3[0]},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	require.Equal(t, 3, fastckks.MaintainedLimbCount(&ckksParams, ckksParams.MaxLevel()))
	return Parameters{
		BootstrappingParameters: ckksParams,
		Mod1ParametersLiteral: mod1.ParametersLiteral{
			LogScale:        30,
			LogMessageRatio: 2,
		},
	}, ckksParams
}

func newFastModUpInput(params ckks.Parameters, values []int64, scale rlwe.Scale, withBacking bool) *rlwe.Ciphertext {
	level := 0
	if withBacking {
		level = params.MaxLevel()
	}
	ct := newFastModUpInputAtLevel(params, values, scale, level)
	if withBacking {
		ct.Resize(ct.Degree(), 0)
	}
	return ct
}

func newFastModUpInputAtLevel(params ckks.Parameters, values []int64, scale rlwe.Scale, level int) *rlwe.Ciphertext {
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
		if level > 0 {
			for limb := 1; limb <= level; limb++ {
				for j := range ct.Value[component].Coeffs[limb] {
					ct.Value[component].Coeffs[limb][j] = uint64(0x9e3779b9) + uint64(limb*params.N()+j+component)
				}
			}
		}
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

// canonicalModUpBasisReference independently implements the accepted odd-q0
// centered convention, including the positive midpoint r=q0>>1.
func canonicalModUpBasisReference(params Parameters, ct *rlwe.Ciphertext) *rlwe.Ciphertext {
	ringQ := params.BootstrappingParameters.RingQ()
	for component := range ct.Value {
		ringQ.AtLevel(0).INTT(ct.Value[component], ct.Value[component])
	}
	ct.Resize(ct.Degree(), params.BootstrappingParameters.MaxLevel())
	q0 := new(big.Int).SetUint64(ringQ.SubRings[0].Modulus)
	half := ringQ.SubRings[0].Modulus / 2
	for component := range ct.Value {
		for coefficient, residue := range ct.Value[component].Coeffs[0] {
			canonical := new(big.Int).SetUint64(residue)
			if residue > half {
				canonical.Sub(canonical, q0)
			}
			for logicalIndex := 1; logicalIndex <= params.BootstrappingParameters.MaxLevel(); logicalIndex++ {
				modulus := new(big.Int).SetUint64(ringQ.SubRings[logicalIndex].Modulus)
				ct.Value[component].Coeffs[logicalIndex][coefficient] = new(big.Int).Mod(new(big.Int).Set(canonical), modulus).Uint64()
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

func standardModUpTraceReference(params Parameters, ct *rlwe.Ciphertext, logN int) {
	ringQ := params.BootstrappingParameters.RingQ().AtLevel(ct.Level())
	gap := 1 << (params.BootstrappingParameters.LogN() - logN - 1)
	if logN == 0 {
		gap <<= 1
	}
	if gap <= 1 {
		return
	}
	nInv := new(big.Int).SetUint64(uint64(gap))
	if nInv.ModInverse(nInv, ringQ.ModulusAtLevel[ct.Level()]) == nil {
		panic("ModUp trace reference inverse does not exist")
	}
	for component := range ct.Value {
		ringQ.MulScalarBigint(ct.Value[component], nInv, ct.Value[component])
	}
	apply := func(galEl uint64) {
		index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galEl)
		if err != nil {
			panic(err)
		}
		for component := range ct.Value {
			tmp := ringQ.NewPoly()
			ringQ.AutomorphismNTTWithIndex(ct.Value[component], index, tmp)
			ringQ.Add(ct.Value[component], tmp, ct.Value[component])
		}
	}
	for i := logN; i < params.BootstrappingParameters.LogN()-1; i++ {
		apply(params.BootstrappingParameters.GaloisElement(1 << i))
	}
	if logN == 0 {
		apply(ringQ.NthRoot() - 1)
	}
}

func standardCompleteModUpReference(params Parameters, ct *rlwe.Ciphertext) {
	standardModUpBasisReference(params, ct)
	standardModUpTraceReference(params, ct, params.CoeffsToSlotsParameters.LogSlots)
	ringQ := params.BootstrappingParameters.RingQ()
	for component := range ct.Value {
		ringQ.MForm(ct.Value[component], ct.Value[component])
	}
	ct.IsMontgomery = true
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

func requireFastModUpMaintainedMatches(t *testing.T, params Parameters, want, got *rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, want.Level(), got.Level())
	require.Equal(t, want.Degree(), got.Degree())
	require.True(t, want.Scale.Equal(got.Scale))
	require.Equal(t, want.IsNTT, got.IsNTT)
	require.Equal(t, want.IsMontgomery, got.IsMontgomery)
	require.Equal(t, want.IsBatched, got.IsBatched)
	require.Equal(t, want.IsBitReversed, got.IsBitReversed)
	require.Equal(t, want.LogDimensions, got.LogDimensions)
	maintained := fastckks.MaintainedLimbCount(&params.BootstrappingParameters, got.Level())
	for component := range got.Value {
		for logicalIndex := 0; logicalIndex < maintained; logicalIndex++ {
			require.Equal(t, want.Value[component].Coeffs[logicalIndex], got.Value[component].Coeffs[logicalIndex],
				"component=%d q-index=%d", component, logicalIndex)
		}
	}
}

func unfusedPrivateFModUpBasisReference(t *testing.T, eval *FastEvaluator, source *rlwe.Ciphertext) *rlwe.Ciphertext {
	t.Helper()
	params := eval.Parameters.BootstrappingParameters
	private, err := fastckks.ImportLevel0(params, source, 3, fastckks.FastCiphertextDomain{IsNTT: source.IsNTT})
	require.NoError(t, err)
	private, err = fastckks.FastStorageModUpLevel0(private, params.MaxLevel())
	require.NoError(t, err)
	compact, err := private.ExportToCompactLogical(fastckks.FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	if scale := (eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()) / compact.Scale.Float64(); scale > 1 {
		scalar := uint64(math.Round(scale))
		err = eval.FastCKKS.MulIntegerMaintained(compact, new(big.Int).SetUint64(scalar), compact)
		require.NoError(t, err)
		compact.Scale = compact.Scale.Mul(rlwe.NewScale(scale))
	}
	return compact
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
	input.IsBatched = false
	input.IsBitReversed = true
	input.LogDimensions = ring.Dimensions{Rows: 1, Cols: 2}
	canonicalWant := canonicalModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<20), false))
	canonicalWant.IsBatched = input.IsBatched
	canonicalWant.IsBitReversed = input.IsBitReversed
	canonicalWant.LogDimensions = input.LogDimensions
	want := unfusedPrivateFModUpBasisReference(t, fastEval, input)
	requireFastModUpBasisMatches(t, canonicalWant, want)
	got, err := fastEval.modUpBasis(input)
	require.NoError(t, err)
	requireFastModUpBasisMatches(t, want, got)
	require.Equal(t, params.BootstrappingParameters.MaxLevel(), got.Level())
	require.True(t, got.IsNTT)
	require.False(t, got.IsMontgomery)
	require.Equal(t, input.IsBatched, got.IsBatched)
	require.Equal(t, input.IsBitReversed, got.IsBitReversed)
	require.Equal(t, input.LogDimensions, got.LogDimensions)
}

func TestFastModUpBasisQ012CanonicalResiduesAndMetadata(t *testing.T) {
	params, ckksParams := fastModUpQ012Parameters(t)
	q0 := ckksParams.RingQ().SubRings[0].Modulus
	values := []int64{int64(q0 >> 1), int64(q0>>1) - 1, -int64(q0 >> 1), 12345, -67890, 0}
	input := newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<20), true)
	input.IsBatched = false
	input.IsBitReversed = true
	input.LogDimensions = ring.Dimensions{Rows: 2, Cols: 3}
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	canonicalWant := canonicalModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<20), false))
	canonicalWant.IsBatched = input.IsBatched
	canonicalWant.IsBitReversed = input.IsBitReversed
	canonicalWant.LogDimensions = input.LogDimensions
	want := unfusedPrivateFModUpBasisReference(t, fastEval, input)
	requireFastModUpMaintainedMatches(t, params, canonicalWant, want)
	got, err := fastEval.modUpBasis(input)
	require.NoError(t, err)
	requireFastModUpMaintainedMatches(t, params, want, got)
	require.True(t, got.IsNTT)
	require.False(t, got.IsMontgomery)
	for component := range got.Value {
		for logicalIndex := fastckks.MaintainedLimbCount(&ckksParams, got.Level()); logicalIndex <= got.Level(); logicalIndex++ {
			require.Nil(t, got.Value[component].Coeffs[logicalIndex], "dormant q-index %d", logicalIndex)
		}
	}
}

func TestFastModUpBasisMatchesHistoricalAwayFromMidpoint(t *testing.T) {
	params, ckksParams := fastModUpParameters(t, 4)
	q0 := ckksParams.RingQ().SubRings[0].Modulus
	values := []int64{int64(q0>>1) - 1, int64(q0>>1) + 1, -int64(q0 >> 1), 12345, -67890, 0}
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	got, err := fastEval.modUpBasis(newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<20), false))
	require.NoError(t, err)
	want := standardModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<20), false))
	requireFastModUpBasisMatches(t, want, got)
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
	q0Rows := make([][]uint64, 2)
	for component := range input.Value {
		for j := range input.Value[component].Coeffs[0] {
			value := values[(j+component)%len(values)]
			input.Value[component].Coeffs[0][j] = new(big.Int).Mod(big.NewInt(value), new(big.Int).SetUint64(q0)).Uint64()
		}
		ckksParams.RingQ().SubRings[0].NTT(input.Value[component].Coeffs[0], input.Value[component].Coeffs[0])
		q0Rows[component] = append([]uint64(nil), input.Value[component].Coeffs[0]...)
		for limb := 2; limb <= ckksParams.MaxLevel(); limb++ {
			for j := range input.Value[component].Coeffs[limb] {
				input.Value[component].Coeffs[limb][j] = uint64(0x9e3779b9) + uint64(limb*ckksParams.N()+j+component)
			}
		}
	}
	input.Resize(input.Degree(), 0)
	want := canonicalModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<40), false))
	got, err := fastEval.modUpBasis(input)
	require.NoError(t, err)
	requireFastModUpBasisMatches(t, want, got)
	require.Equal(t, rlwe.NewScale(1<<40), got.Scale)
	for component := range input.Value {
		require.Equal(t, q0Rows[component], input.Value[component].Coeffs[0], "input q0 NTT row must not be mutated")
	}
	for component := range got.Value {
		for limb := 2; limb <= ckksParams.MaxLevel(); limb++ {
			require.Nil(t, got.Value[component].Coeffs[limb], "dormant output row %d must not be materialized", limb)
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
	want := canonicalModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<40), false))
	requireFastModUpBasisMatches(t, want, got)
}

func TestFastModUpCompleteBoundary(t *testing.T) {
	values := []int64{1, -2, 12345, -67890}
	for _, profile := range []struct {
		name string
		q012 bool
	}{
		{name: "q01"},
		{name: "q012", q012: true},
	} {
		t.Run(profile.name, func(t *testing.T) {
			var params Parameters
			var ckksParams ckks.Parameters
			if profile.q012 {
				params, ckksParams = fastModUpQ012Parameters(t)
			} else {
				params, ckksParams = fastModUpParameters(t, 4)
			}
			maintained := fastckks.MaintainedLimbCount(&ckksParams, ckksParams.MaxLevel())
			for _, logSlots := range []int{3, 1} {
				t.Run(fmt.Sprintf("logSlots-%d", logSlots), func(t *testing.T) {
					params.CoeffsToSlotsParameters.LogSlots = logSlots
					fastEval, err := NewFastEvaluator(params)
					require.NoError(t, err)
					input := newFastModUpInputAtLevel(ckksParams, values, rlwe.NewScale(1<<20), ckksParams.MaxLevel())
					input.Resize(input.Degree(), 0)

					want := canonicalModUpBasisReference(params, newFastModUpInput(ckksParams, values, rlwe.NewScale(1<<20), false))
					standardModUpTraceReference(params, want, logSlots)
					for component := range want.Value {
						for limb := 0; limb < maintained; limb++ {
							ckksParams.RingQ().SubRings[limb].MForm(want.Value[component].Coeffs[limb], want.Value[component].Coeffs[limb])
						}
					}
					want.IsMontgomery = true

					got, err := fastEval.ModUp(input)
					require.NoError(t, err)
					require.Equal(t, want.Level(), got.Level())
					require.True(t, want.Scale.Equal(got.Scale))
					require.True(t, got.IsNTT)
					require.True(t, got.IsMontgomery)
					for component := range got.Value {
						for limb := 0; limb < maintained; limb++ {
							require.Equal(t, want.Value[component].Coeffs[limb], got.Value[component].Coeffs[limb], "component %d limb %d", component, limb)
						}
						for limb := maintained; limb <= ckksParams.MaxLevel(); limb++ {
							require.Nil(t, got.Value[component].Coeffs[limb], "Trace/ModUp output must keep dormant row %d unmaterialized", limb)
						}
					}
				})
			}
		})
	}
}

func TestFastStorageNormalizedTraceMatchesLogicalTrace(t *testing.T) {
	values := []int64{1, -2, 12345, -67890}
	for _, profile := range []struct {
		name string
		new  func(testing.TB) (Parameters, ckks.Parameters)
	}{
		{name: "q01", new: func(tb testing.TB) (Parameters, ckks.Parameters) { return fastModUpParameters(tb, 4) }},
		{name: "q012", new: fastModUpQ012Parameters},
	} {
		t.Run(profile.name, func(t *testing.T) {
			params, ckksParams := profile.new(t)
			maintained := fastckks.MaintainedLimbCount(&ckksParams, ckksParams.MaxLevel())
			for _, logSlots := range []int{3, 2, 1, 0} {
				for _, scale := range []int64{1 << 20, 1 << 28} {
					t.Run(fmt.Sprintf("logSlots-%d-scale-%d", logSlots, scale), func(t *testing.T) {
						input := newFastModUpInput(ckksParams, values, rlwe.NewScale(scale), false)
						fastEvaluator := fastckks.NewEvaluator(ckksParams)
						reference, err := fastckks.FusedLevel0ModUpToCompactLogical(ckksParams, input, ckksParams.MaxLevel(), fastckks.FastCiphertextDomain{IsNTT: true})
						require.NoError(t, err)
						private, err := fastckks.ImportLevel0(ckksParams, input, 3, fastckks.FastCiphertextDomain{IsNTT: true})
						require.NoError(t, err)
						private, err = fastckks.FastStorageModUpLevel0(private, ckksParams.MaxLevel())
						require.NoError(t, err)

						mod1Parameters := mod1.Parameters{LogDefaultScale: params.Mod1ParametersLiteral.LogScale, LogMessageRatio: params.Mod1ParametersLiteral.LogMessageRatio}
						if alignment := (mod1Parameters.ScalingFactor().Float64() / mod1Parameters.MessageRatio()) / private.Scale().Float64(); alignment > 1 {
							scalar := uint64(math.Round(alignment))
							private, err = fastckks.FastStorageMulInteger(private, new(big.Int).SetUint64(scalar))
							require.NoError(t, err)
							require.NoError(t, private.SetScale(private.Scale().Mul(rlwe.NewScale(alignment))))
							require.NoError(t, fastEvaluator.MulIntegerMaintained(reference, new(big.Int).SetUint64(scalar), reference))
							reference.Scale = reference.Scale.Mul(rlwe.NewScale(alignment))
						}

						require.NoError(t, fastEvaluator.Trace(reference, logSlots, reference))
						private, err = fastckks.FastStorageTraceNormalized(private, logSlots)
						require.NoError(t, err)
						got, err := private.ExportToCompactLogical(fastckks.FastCiphertextDomain{IsNTT: true})
						require.NoError(t, err)

						require.Equal(t, reference.Level(), got.Level())
						require.True(t, reference.Scale.Equal(got.Scale))
						require.Equal(t, reference.IsNTT, got.IsNTT)
						require.Equal(t, reference.IsMontgomery, got.IsMontgomery)
						require.Equal(t, reference.MetaData.IsBatched, got.MetaData.IsBatched)
						require.Equal(t, reference.MetaData.IsBitReversed, got.MetaData.IsBitReversed)
						require.Equal(t, reference.MetaData.LogDimensions, got.MetaData.LogDimensions)
						for component := range reference.Value {
							for limb := 0; limb < maintained; limb++ {
								require.Equal(t, reference.Value[component].Coeffs[limb], got.Value[component].Coeffs[limb], "component=%d limb=%d", component, limb)
							}
							for limb := maintained; limb <= ckksParams.MaxLevel(); limb++ {
								require.Nil(t, got.Value[component].Coeffs[limb], "dormant q row %d must remain unmaterialized", limb)
							}
						}
					})
				}
			}
		})
	}
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
	ct.IsNTT = true
	ct.IsMontgomery = false
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

func BenchmarkFastModUpLogN13(b *testing.B) {
	benchmarkFastModUp(b, 13, true)
}

func BenchmarkStandardModUpReferenceLogN13(b *testing.B) {
	benchmarkFastModUp(b, 13, false)
}

func BenchmarkFastModUpLogN16(b *testing.B) {
	benchmarkFastModUp(b, 16, true)
}

func BenchmarkStandardModUpReferenceLogN16(b *testing.B) {
	benchmarkFastModUp(b, 16, false)
}

func benchmarkFastModUp(b *testing.B, logN int, fastPath bool) {
	params, ckksParams := fastModUpParameters(b, logN)
	params.CoeffsToSlotsParameters.LogSlots = 1
	values := []int64{1 << 40, -(1 << 40) + 17, 1 << 39, -123456789}
	scale := rlwe.NewScale(1 << 40)
	template := newFastModUpInput(ckksParams, values, scale, true)
	input := newFastModUpInput(ckksParams, values, scale, true)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(b, err)
	if fastPath {
		if _, err := fastEval.ModUp(input); err != nil {
			b.Fatal(err)
		}
		resetFastModUpInput(input, template, scale)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		resetFastModUpInput(input, template, scale)
		b.StartTimer()
		if fastPath {
			if _, err := fastEval.ModUp(input); err != nil {
				b.Fatal(err)
			}
		} else {
			standardCompleteModUpReference(params, input)
		}
	}
}
