package bootstrapping

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func fastScaleDownParameters(t testing.TB, logN int) (Parameters, ckks.Parameters) {
	t.Helper()
	logQ := []int{55, 39, 50}
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
			LogMessageRatio: 2,
		},
	}, params
}

func fastScaleDownActualLogN16Parameters(t testing.TB) (Parameters, ckks.Parameters) {
	t.Helper()
	logQ := []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            16,
		LogQ:            logQ,
		LogDefaultScale: 45,
	})
	require.NoError(t, err)
	require.Equal(t, 56, new(big.Int).SetUint64(params.Q()[0]).BitLen())
	require.Equal(t, 39, new(big.Int).SetUint64(params.Q()[1]).BitLen())
	return Parameters{
		BootstrappingParameters: params,
		Mod1ParametersLiteral: mod1.ParametersLiteral{
			LogMessageRatio: 10,
		},
	}, params
}

func scaleQ0Div(params ckks.Parameters, divisor int64) rlwe.Scale {
	value := new(big.Float).SetPrec(256).SetUint64(params.Q()[0])
	value.Quo(value, new(big.Float).SetInt64(divisor))
	return rlwe.NewScale(value)
}

func newScaleDownCiphertext(params ckks.Parameters, level int, scale rlwe.Scale, values []int64) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, level)
	ct.IsNTT = true
	ct.Scale = scale
	for d := range ct.Value {
		for limb := 0; limb <= level; limb++ {
			modulus := new(big.Int).SetUint64(params.RingQ().SubRings[limb].Modulus)
			for i := range ct.Value[d].Coeffs[limb] {
				value := values[(i+d)%len(values)]
				if d == 1 {
					value = -value
				}
				ct.Value[d].Coeffs[limb][i] = new(big.Int).Mod(big.NewInt(value), modulus).Uint64()
			}
		}
		params.RingQ().AtLevel(level).NTT(ct.Value[d], ct.Value[d])
	}
	return ct
}

func poisonScaleDownDormant(ct *rlwe.Ciphertext) {
	for d := range ct.Value {
		for limb := 2; limb <= ct.Level(); limb++ {
			for i := range ct.Value[d].Coeffs[limb] {
				ct.Value[d].Coeffs[limb][i] = ^uint64(0) - uint64(i+limb+d)
			}
		}
	}
}

func standardScaleDownReference(params Parameters, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, *rlwe.Scale, error) {
	eval := ckks.NewEvaluator(params.BootstrappingParameters, nil)
	r := params.BootstrappingParameters.RingQ()
	msgRatio := float64(uint(1 << params.Mod1ParametersLiteral.LogMessageRatio))
	for ct.Level() != 0 && checkMessageRatio(ct, msgRatio, r) {
		ct.Resize(ct.Degree(), ct.Level()-1)
	}
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[ct.Level()]).Div(ct.Scale)
	targetMessageRatio := rlwe.NewScale(msgRatio)
	scaleUp := currentMessageRatio.Div(targetMessageRatio)
	if scaleUp.Cmp(rlwe.NewScale(0.5)) == -1 {
		return nil, nil, nil
	}
	scaleUpBigint := scaleUp.BigInt()
	if err := eval.Mul(ct, scaleUpBigint, ct); err != nil {
		return nil, nil, err
	}
	ct.Scale = ct.Scale.Mul(rlwe.NewScale(scaleUpBigint))
	targetScale := new(big.Float).SetPrec(256).SetInt(r.ModulusAtLevel[0])
	targetScale.Quo(targetScale, new(big.Float).SetFloat64(msgRatio))
	if ct.Level() != 0 {
		if err := eval.RescaleTo(ct, rlwe.NewScale(targetScale), ct); err != nil {
			return nil, nil, err
		}
	}
	errScale := ct.Scale.Div(rlwe.NewScale(targetScale))
	return ct, &errScale, nil
}

func requireScaleDownMatches(t *testing.T, want, got *rlwe.Ciphertext, wantErr, gotErr *rlwe.Scale) {
	t.Helper()
	require.Equal(t, 0, got.Level())
	require.Equal(t, 0, want.Level())
	require.True(t, want.Scale.Equal(got.Scale))
	require.True(t, (*wantErr).Equal(*gotErr))
	for d := range want.Value {
		require.Equal(t, want.Value[d].Coeffs[0], got.Value[d].Coeffs[0], "component %d r0", d)
	}
}

func TestFastScaleDownLevelZero(t *testing.T) {
	params, ckksParams := fastScaleDownParameters(t, 4)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	source := newScaleDownCiphertext(ckksParams, 0, scaleQ0Div(ckksParams, 8), []int64{1 << 40, -(1 << 40) + 17, 12345})
	want, wantErr, err := standardScaleDownReference(params, source.CopyNew())
	require.NoError(t, err)
	got, gotErr, err := fastEval.ScaleDown(source)
	require.NoError(t, err)
	requireScaleDownMatches(t, want, got, wantErr, gotErr)
}

func TestFastScaleDownCheapDropsToLevelZeroAndPoison(t *testing.T) {
	params, ckksParams := fastScaleDownParameters(t, 4)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	source := newScaleDownCiphertext(ckksParams, 2, scaleQ0Div(ckksParams, 4), []int64{1 << 40, -(1 << 40) + 31, 777777})
	want, wantErr, err := standardScaleDownReference(params, source.CopyNew())
	require.NoError(t, err)
	poisoned := source.CopyNew()
	poisonScaleDownDormant(poisoned)
	got, gotErr, err := fastEval.ScaleDown(poisoned)
	require.NoError(t, err)
	requireScaleDownMatches(t, want, got, wantErr, gotErr)
	require.Equal(t, 0, got.Level())
}

func TestFastScaleDownRescalesFromLevelOneAndMatchesStandard(t *testing.T) {
	params, ckksParams := fastScaleDownParameters(t, 4)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	source := newScaleDownCiphertext(ckksParams, 2, scaleQ0Div(ckksParams, 2), []int64{1 << 40, -(1 << 40) + 43, 987654321})
	want, wantErr, err := standardScaleDownReference(params, source.CopyNew())
	require.NoError(t, err)
	poisoned := source.CopyNew()
	poisonScaleDownDormant(poisoned)
	got, gotErr, err := fastEval.ScaleDown(poisoned)
	require.NoError(t, err)
	requireScaleDownMatches(t, want, got, wantErr, gotErr)
	require.NotEqual(t, uint64(0), got.Value[0].Coeffs[0][0])
}

func TestFastScaleDownMontgomeryLevelZero(t *testing.T) {
	params, ckksParams := fastScaleDownParameters(t, 4)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	source := newScaleDownCiphertext(ckksParams, 0, scaleQ0Div(ckksParams, 8), []int64{1 << 40, -(1 << 40) + 61, 13579})
	for d := range source.Value {
		ckksParams.RingQ().SubRings[0].MForm(source.Value[d].Coeffs[0], source.Value[d].Coeffs[0])
	}
	source.IsMontgomery = true
	want, wantErr, err := standardScaleDownReference(params, source.CopyNew())
	require.NoError(t, err)
	got, gotErr, err := fastEval.ScaleDown(source)
	require.NoError(t, err)
	requireScaleDownMatches(t, want, got, wantErr, gotErr)
}

func TestFastScaleDownActualLogN16Profile(t *testing.T) {
	params, ckksParams := fastScaleDownActualLogN16Parameters(t)
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	inputLevel := ckksParams.MaxLevel()
	source := newScaleDownCiphertext(ckksParams, inputLevel, scaleQ0Div(ckksParams, 2), []int64{1 << 40, -(1 << 40) + 17, 12345})
	want, wantErr, err := standardScaleDownReference(params, source.CopyNew())
	require.NoError(t, err)
	got, gotErr, err := fastEval.ScaleDown(source)
	require.NoError(t, err)
	requireScaleDownMatches(t, want, got, wantErr, gotErr)
}

func BenchmarkFastScaleDown(b *testing.B) {
	for _, logN := range []int{13, 16} {
		params, ckksParams := fastScaleDownParameters(b, logN)
		inputLevel := ckksParams.MaxLevel()
		inputScale := scaleQ0Div(ckksParams, 2)
		values := []int64{1 << 40, -(1 << 40) + 17, 1 << 39, -123456789}
		fastEval, err := NewFastEvaluator(params)
		require.NoError(b, err)
		standardEval := ckks.NewEvaluator(ckksParams, nil)
		fastInput := newScaleDownCiphertext(ckksParams, inputLevel, inputScale, values)
		standardInput := fastInput.CopyNew()

		b.Run("LogN"+big.NewInt(int64(logN)).String()+"/Fast", func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				b.StopTimer()
				fastInput.Resize(1, inputLevel)
				fastInput.Scale = inputScale
				b.StartTimer()
				if _, _, err := fastEval.ScaleDown(fastInput); err != nil {
					b.Fatal(err)
				}
			}
		})
		b.Run("LogN"+big.NewInt(int64(logN)).String()+"/Standard", func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				b.StopTimer()
				standardInput.Resize(1, inputLevel)
				standardInput.Scale = inputScale
				b.StartTimer()
				if _, _, err := standardScaleDownWithEvaluator(params, standardEval, standardInput); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

func standardScaleDownWithEvaluator(params Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, *rlwe.Scale, error) {
	r := params.BootstrappingParameters.RingQ()
	msgRatio := float64(uint(1 << params.Mod1ParametersLiteral.LogMessageRatio))
	for ct.Level() != 0 && checkMessageRatio(ct, msgRatio, r) {
		ct.Resize(ct.Degree(), ct.Level()-1)
	}
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[ct.Level()]).Div(ct.Scale)
	scaleUp := currentMessageRatio.Div(rlwe.NewScale(msgRatio))
	scaleUpBigint := scaleUp.BigInt()
	if err := eval.Mul(ct, scaleUpBigint, ct); err != nil {
		return nil, nil, err
	}
	ct.Scale = ct.Scale.Mul(rlwe.NewScale(scaleUpBigint))
	targetScale := new(big.Float).SetPrec(256).SetInt(r.ModulusAtLevel[0])
	targetScale.Quo(targetScale, new(big.Float).SetFloat64(msgRatio))
	if ct.Level() != 0 {
		if err := eval.RescaleTo(ct, rlwe.NewScale(targetScale), ct); err != nil {
			return nil, nil, err
		}
	}
	errScale := ct.Scale.Div(rlwe.NewScale(targetScale))
	return ct, &errScale, nil
}
