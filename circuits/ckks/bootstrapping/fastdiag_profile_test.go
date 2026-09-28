package bootstrapping

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func fastDiagP93Parameters(t testing.TB) (Parameters, ckks.Parameters) {
	t.Helper()
	residual, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39}, LogDefaultScale: 45, Xs: ring.Ternary{H: 192},
	})
	require.NoError(t, err)
	logN, logSlots, zero, evalModScale, degree, doubleAngle := 13, 12, 0, 60, 30, 3
	k, logMessageRatio := 16, 10
	params, err := NewParametersFromLiteral(residual, ParametersLiteral{
		LogN: &logN, LogSlots: &logSlots, LogP: []int{61, 61, 61, 61, 61},
		SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{39}, {39}, {39}},
		CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{56}, {56}, {56}, {56}},
		EvalModLogScale: &evalModScale, Mod1Degree: &degree, DoubleAngle: &doubleAngle,
		K: &k, LogMessageRatio: &logMessageRatio, Mod1InvDegree: &zero,
		EphemeralSecretWeight: &zero,
	})
	require.NoError(t, err)
	params.CircuitOrder = ModUpThenEncode
	params.ResidualParameters = residual
	return params, residual
}

func fastDiagP93Fixture(t testing.TB) (Parameters, ckks.Parameters, []complex128, *rlwe.Ciphertext) {
	t.Helper()
	params, residual := fastDiagP93Parameters(t)
	values := make([]complex128, 1<<params.CoeffsToSlotsParameters.LogSlots)
	for j := range values {
		values[j] = complex(float64(j%17-8)/256, float64((3*j)%13-6)/512)
	}
	input := fastBootstrapEncodedCiphertext(t, residual, params.CoeffsToSlotsParameters.LogSlots, values)
	return params, residual, values, input
}

func BenchmarkFastDiagP93Q55Count1(b *testing.B) {
	params, residual, _, input := fastDiagP93Fixture(b)
	eval, err := NewFastEvaluator(params)
	require.NoError(b, err)
	if _, err = eval.Bootstrap(input.CopyNew()); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err = eval.Bootstrap(input.CopyNew()); err != nil {
			b.Fatal(err)
		}
	}
	b.StopTimer()
	b.Logf("profile=p93-q55 LogN=%d LogSlots=%d batch=1 residualQ=%v", params.BootstrappingParameters.LogN(), params.CoeffsToSlotsParameters.LogSlots, residual.Q())
}
