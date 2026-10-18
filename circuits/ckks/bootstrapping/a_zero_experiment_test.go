package bootstrapping

import (
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

func requireBootstrapAZero(t *testing.T, evk *rlwe.EvaluationKey) {
	t.Helper()
	if evk == nil {
		return
	}
	require.Equal(t, 1, evk.Degree())
	for i := range evk.Value {
		for j := range evk.Value[i] {
			for level := range evk.Value[i][j][1].Q.Coeffs {
				for _, coefficient := range evk.Value[i][j][1].Q.Coeffs[level] {
					require.Zero(t, coefficient)
				}
			}
			for level := range evk.Value[i][j][1].P.Coeffs {
				for _, coefficient := range evk.Value[i][j][1].P.Coeffs[level] {
					require.Zero(t, coefficient)
				}
			}
		}
	}
}

func TestAZeroBootstrapEvaluationKeySet(t *testing.T) {
	// This constructs the real single-party Bootstrap key set at a small ring
	// degree. No Bootstrap evaluation is run in this key-generation test.
	schemeParams, err := ckks.NewParametersFromLiteral(testPrec45)
	require.NoError(t, err)
	btpLit := ParametersLiteral{
		LogN:                  utils.Pointy(10),
		LogSlots:              utils.Pointy(1),
		EphemeralSecretWeight: utils.Pointy(32),
	}
	btpParams, err := NewParametersFromLiteral(schemeParams, btpLit)
	require.NoError(t, err)
	require.Equal(t, 32, btpParams.EphemeralSecretWeight)
	require.Equal(t, ring.Standard, btpParams.ResidualParameters.RingType())

	sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()
	btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
	require.NoError(t, err)

	requireBootstrapAZero(t, btpKeys.EvkN1ToN2)
	requireBootstrapAZero(t, btpKeys.EvkN2ToN1)
	requireBootstrapAZero(t, btpKeys.EvkRealToCmplx)
	requireBootstrapAZero(t, btpKeys.EvkCmplxToReal)
	require.NotNil(t, btpKeys.EvkDenseToSparse)
	require.NotNil(t, btpKeys.EvkSparseToDense)
	requireBootstrapAZero(t, btpKeys.EvkDenseToSparse)
	requireBootstrapAZero(t, btpKeys.EvkSparseToDense)

	rlk, err := btpKeys.GetRelinearizationKey()
	require.NoError(t, err)
	requireBootstrapAZero(t, &rlk.EvaluationKey)
	for _, galEl := range btpKeys.GetGaloisKeysList() {
		gk, err := btpKeys.GetGaloisKey(galEl)
		require.NoError(t, err)
		requireBootstrapAZero(t, &gk.EvaluationKey)
	}
}
