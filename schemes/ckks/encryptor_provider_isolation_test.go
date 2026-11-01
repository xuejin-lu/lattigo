package ckks

import (
	"math/cmplx"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// Expose only the RLWE provider interface, without Fast's CKKS capability.
// This exercises native CKKS encryption in this fork, not a formal Standard baseline.
type nativeEncryptionParameters struct{ rlwe.ParameterProvider }

func TestCKKSPublicKeyEncryptionWithoutFastCapability(t *testing.T) {
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: 8, LogQ: []int{55, 45}, LogDefaultScale: 40,
	})
	require.NoError(t, err)
	nativeParams := nativeEncryptionParameters{params}
	kgen := rlwe.NewKeyGenerator(nativeParams)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	pt := NewPlaintext(params, 1)
	encoder := NewEncoder(params)
	values := make([]complex128, params.MaxSlots())
	for i := range values {
		values[i] = complex(float64(i%11)/32, -float64(i%7)/64)
	}
	require.NoError(t, encoder.Encode(values, pt))
	ct, err := rlwe.NewEncryptor(nativeParams, pk).EncryptNew(pt)
	require.NoError(t, err)
	nonzeroC1 := false
	for _, row := range ct.Value[1].Coeffs {
		for _, coefficient := range row {
			nonzeroC1 = nonzeroC1 || coefficient != 0
		}
	}
	require.True(t, nonzeroC1, "a provider without Fast capability must retain native PK encryption")
	require.False(t, pt.Value.Equal(&ct.Value[0]), "native PK encryption must mask c0")
	require.True(t, pt.MetaData.Equal(ct.MetaData))
	decoded := make([]complex128, len(values))
	require.NoError(t, encoder.Decode(rlwe.NewDecryptor(nativeParams, sk).DecryptNew(ct), decoded))
	for i := range values {
		require.Less(t, cmplx.Abs(decoded[i]-values[i]), 1e-7)
	}
}
