package bootstrapping

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

func TestFastBootstrapNoiseConfigurationForwardsWithoutRunningBootstrap(t *testing.T) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39, 40, 39, 40, 60, 60}, LogP: []int{61},
		LogDefaultScale: 45, RingType: ring.Standard,
	})
	require.NoError(t, err)
	fastEval := &FastEvaluator{FastCKKS: fastckks.NewEvaluator(params)}
	eval := &Evaluator{fast: fastEval}

	initial, ok := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot().(fastckks.FixedSigmaBootstrapNoiseSnapshot)
	require.True(t, ok)
	require.False(t, initial.Enabled)
	require.NoError(t, eval.EnableFixedSigmaBootstrapNoiseDiagnostic(4.004280012097432, 20261027))
	configured, ok := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot().(fastckks.FixedSigmaBootstrapNoiseSnapshot)
	require.True(t, ok)
	require.True(t, configured.Enabled)
	require.Equal(t, 4.004280012097432, configured.Sigma)
	require.Equal(t, int64(20261027), configured.RootSeed)
	eval.DisableFixedSigmaBootstrapNoiseDiagnostic()
	disabled, ok := eval.FixedSigmaBootstrapNoiseDiagnosticSnapshot().(fastckks.FixedSigmaBootstrapNoiseSnapshot)
	require.True(t, ok)
	require.False(t, disabled.Enabled)
}
