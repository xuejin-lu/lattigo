package bootstrapping

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

func fastBootstrapParameters(t testing.TB, logSlots int, factorTwo bool) (Parameters, ckks.Parameters) {
	t.Helper()
	logN := 4
	return fastBootstrapParametersAt(t, logN, logSlots, factorTwo)
}

func fastBootstrapParametersAt(t testing.TB, logN, logSlots int, factorTwo bool) (Parameters, ckks.Parameters) {
	t.Helper()
	residualLogN := logN
	if factorTwo {
		residualLogN = logN - 1
	}
	residualLiteral := ckks.ParametersLiteral{LogN: residualLogN, LogQ: []int{55, 39}, LogDefaultScale: 30}
	if logN == 16 {
		residualLiteral.LogQ = []int{54, 38}
	}
	if factorTwo {
		residualLiteral.LogNthRoot = logN + 1
	}
	residual, err := ckks.NewParametersFromLiteral(residualLiteral)
	require.NoError(t, err)
	zero := 0
	factorization := [][]int{{39}}
	if logSlots > 1 {
		factorization = [][]int{{39}, {39}}
	}
	literal := ParametersLiteral{
		LogN:     utils.Pointy(logN),
		LogSlots: utils.Pointy(logSlots),
		SlotsToCoeffsFactorizationDepthAndLogScales: factorization,
		CoeffsToSlotsFactorizationDepthAndLogScales: factorization,
		EvalModLogScale:       utils.Pointy(45),
		Mod1Degree:            utils.Pointy(30),
		DoubleAngle:           utils.Pointy(2),
		K:                     utils.Pointy(4),
		LogMessageRatio:       utils.Pointy(2),
		Mod1InvDegree:         &zero,
		EphemeralSecretWeight: &zero,
	}
	params, err := NewParametersFromLiteral(residual, literal)
	require.NoError(t, err)
	return params, residual
}

func fastBootstrapPlainCiphertext(t testing.TB, params ckks.Parameters, level, logSlots int) *rlwe.Ciphertext {
	t.Helper()
	ct := ckks.NewCiphertext(params, 1, level)
	ct.IsNTT = true
	ct.IsMontgomery = false
	ct.IsBatched = true
	ct.Scale = params.DefaultScale()
	ct.LogDimensions = ring.Dimensions{Cols: logSlots}
	for d := range ct.Value {
		for limb := 0; limb <= level; limb++ {
			for i := range ct.Value[d].Coeffs[limb] {
				value := int64((i % 5) - 2)
				if d == 1 {
					value = 0
				}
				ct.Value[d].Coeffs[limb][i] = new(big.Int).Mod(big.NewInt(value), new(big.Int).SetUint64(params.RingQ().SubRings[limb].Modulus)).Uint64()
			}
		}
		params.RingQ().AtLevel(level).NTT(ct.Value[d], ct.Value[d])
	}
	return ct
}

func fastBootstrapEncodedCiphertext(t testing.TB, params ckks.Parameters, logSlots int, values []complex128) *rlwe.Ciphertext {
	return fastBootstrapEncodedCiphertextAtLevel(t, params, 0, logSlots, values)
}

func fastBootstrapEncodedCiphertextAtLevel(t testing.TB, params ckks.Parameters, level, logSlots int, values []complex128) *rlwe.Ciphertext {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	pt := ckks.NewPlaintext(params, level)
	pt.IsNTT = true
	pt.IsMontgomery = false
	pt.LogDimensions = ring.Dimensions{Cols: logSlots}
	require.NoError(t, encoder.Encode(values, pt))
	ct := ckks.NewCiphertext(params, 1, level)
	*ct.MetaData = *pt.MetaData
	ct.Value[0].Copy(pt.Value)
	ct.Value[1].Zero()
	ct.IsNTT = pt.IsNTT
	ct.IsMontgomery = pt.IsMontgomery
	return ct
}

func fastBootstrapDecode(t *testing.T, params ckks.Parameters, ct *rlwe.Ciphertext, values []complex128) {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	pt := ckks.NewPlaintext(params, ct.Level())
	*pt.MetaData = *ct.MetaData
	pt.Value.Copy(ct.Value[0])
	pt.IsNTT = ct.IsNTT
	pt.IsMontgomery = ct.IsMontgomery
	decoded := make([]complex128, len(values))
	require.NoError(t, encoder.Decode(pt, decoded))
	for i := range values {
		require.InDelta(t, real(values[i]), real(decoded[i]), 1e-2, "real slot %d", i)
		require.InDelta(t, imag(values[i]), imag(decoded[i]), 1e-2, "imag slot %d", i)
	}
}

func TestFastBootstrapCircuitInitializationAndInterface(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.Equal(t, 0, eval.MinimumInputLevel())
	require.Equal(t, residual.MaxLevel(), eval.OutputLevel())
	require.NotNil(t, eval)
	require.NoError(t, eval.ensureFastBootstrapCircuit())
	require.NotNil(t, eval.DFTEvaluator)
	require.NotEmpty(t, eval.C2SDFTMatrix.Matrices)
	require.NotEmpty(t, eval.S2CDFTMatrix.Matrices)
	require.Equal(t, mod1.CosDiscrete, eval.Mod1Parameters.Mod1Type)
	require.Nil(t, eval.Mod1Parameters.Mod1InvPoly)
}

func TestFastBootstrapPublicValidation(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	_, err = eval.Bootstrap(nil)
	require.Error(t, err)
	_, err = eval.BootstrapMany(nil)
	require.Error(t, err)

	ct := fastBootstrapPlainCiphertext(t, residual, 0, 2)
	ct.IsMontgomery = true
	_, err = eval.Bootstrap(ct)
	require.Error(t, err)
}

func TestFastBootstrapStageAEndToEndSmoke(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	ct := fastBootstrapPlainCiphertext(t, residual, 0, 2)
	source := ct.CopyNew()
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	require.Equal(t, residual.N(), out.N())
	require.Equal(t, residual.MaxLevel(), out.Level())
	require.True(t, out.IsNTT)
	require.False(t, out.IsMontgomery)
	require.True(t, out.Scale.Equal(residual.DefaultScale()))
	require.Equal(t, source.Value[0].Coeffs[0], ct.Value[0].Coeffs[0])
}

func TestFastBootstrapPreservesNonZeroMessage(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i}
	ct := fastBootstrapEncodedCiphertext(t, residual, 2, values)
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	fastBootstrapDecode(t, residual, out, values)
}

func TestFastBootstrapLevelOneInput(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	ct := fastBootstrapEncodedCiphertextAtLevel(t, residual, 1, 2, []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i})
	source := ct.CopyNew()
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	require.Equal(t, residual.MaxLevel(), out.Level())
	require.Equal(t, source.Value[0].Coeffs[:2], ct.Value[0].Coeffs[:2])
}

func TestFastBootstrapSparseRepackedPath(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 1, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i}
	ct := fastBootstrapEncodedCiphertext(t, residual, 1, values)
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	fastBootstrapDecode(t, residual, out, values)
}

func TestFastBootstrapManyOddCount(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	cts := make([]rlwe.Ciphertext, 3)
	for i := range cts {
		values := []complex128{complex(0.05*float64(i+1), 0.01), complex(-0.04*float64(i+1), -0.02), complex(0.03, 0.05), complex(-0.02, -0.04)}
		cts[i] = *fastBootstrapEncodedCiphertext(t, residual, 2, values)
	}
	outputs, err := eval.BootstrapMany(cts)
	require.NoError(t, err)
	require.Len(t, outputs, 3)
	for _, out := range outputs {
		require.Equal(t, 2, out.LogSlots())
		require.False(t, out.IsMontgomery)
		require.Equal(t, residual.MaxLevel(), out.Level())
	}
}

func TestFastBootstrapMatchesStandardDecodedReference(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i}
	fastInput := fastBootstrapEncodedCiphertext(t, residual, 2, values)
	fastOut, err := fastEval.Bootstrap(fastInput)
	require.NoError(t, err)

	sk := rlwe.NewKeyGenerator(params.BootstrappingParameters).GenSecretKeyNew()
	keys, _, err := params.GenEvaluationKeys(sk)
	require.NoError(t, err)
	standardEval, err := NewEvaluator(params, keys)
	require.NoError(t, err)
	standardInput := fastBootstrapEncodedCiphertext(t, residual, 2, values)
	standardOut, err := standardEval.Bootstrap(standardInput)
	require.NoError(t, err)
	decryptor := rlwe.NewDecryptor(residual, sk)
	standardPlaintext := decryptor.DecryptNew(standardOut)
	standardValues := make([]complex128, len(values))
	require.NoError(t, ckks.NewEncoder(residual).Decode(standardPlaintext, standardValues))
	fastPlaintext := ckks.NewPlaintext(residual, fastOut.Level())
	*fastPlaintext.MetaData = *fastOut.MetaData
	fastPlaintext.Value.Copy(fastOut.Value[0])
	fastPlaintext.IsNTT = fastOut.IsNTT
	fastPlaintext.IsMontgomery = fastOut.IsMontgomery
	fastValues := make([]complex128, len(values))
	require.NoError(t, ckks.NewEncoder(residual).Decode(fastPlaintext, fastValues))
	for i := range values {
		require.InDelta(t, real(standardValues[i]), real(fastValues[i]), 1e-2, "real slot %d", i)
		require.InDelta(t, imag(standardValues[i]), imag(fastValues[i]), 1e-2, "imag slot %d", i)
	}
}

func TestFastBootstrapCoreIgnoresDormantResidues(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.GreaterOrEqual(t, params.BootstrappingParameters.MaxLevel(), 2)
	first := fastBootstrapPlainCiphertext(t, params.BootstrappingParameters, 2, 2)
	q01Scale := new(big.Float).SetUint64(params.BootstrappingParameters.Q()[0])
	q01Scale.Mul(q01Scale, new(big.Float).SetUint64(params.BootstrappingParameters.Q()[1]))
	first.Scale = rlwe.NewScale(q01Scale)
	second := first.CopyNew()
	for d := range second.Value {
		for limb := 2; limb < len(second.Value[d].Coeffs); limb++ {
			for i := range second.Value[d].Coeffs[limb] {
				second.Value[d].Coeffs[limb][i] = uint64(i + 17 + d + limb)
			}
		}
	}
	firstOut, _, err := eval.bootstrapCore(first)
	require.NoError(t, err)
	secondOut, _, err := eval.bootstrapCore(second)
	require.NoError(t, err)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb <= secondOut.Level() && limb < 2; limb++ {
			require.Equal(t, firstOut.Value[d].Coeffs[limb], secondOut.Value[d].Coeffs[limb], "component %d limb %d", d, limb)
		}
	}
}
