package rlwe

import (
	"testing"

	"github.com/stretchr/testify/require"
)

func TestRLWEParameterProviderKeepsNativeEncryptNew(t *testing.T) {
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN:         8,
		LogQ:         []int{55, 45},
		DefaultScale: NewScale(1),
		NTTFlag:      true,
	})
	require.NoError(t, err)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	pt := NewPlaintext(params, 0)
	pt.Value.Coeffs[0][0] = 1

	ct, err := NewEncryptor(params, sk).EncryptNew(pt)
	require.NoError(t, err)
	require.Equal(t, 1, ct.Degree())
	nonzeroC1 := false
	for _, row := range ct.Value[1].Coeffs {
		for _, coefficient := range row {
			nonzeroC1 = nonzeroC1 || coefficient != 0
		}
	}
	require.True(t, nonzeroC1, "generic rlwe.Parameters must retain ordinary RLWE encryption semantics")
}
