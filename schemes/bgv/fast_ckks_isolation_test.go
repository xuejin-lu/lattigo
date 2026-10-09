package bgv

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

func TestBGVSecretEncryptNewRemainsNative(t *testing.T) {
	literal := testInsecure
	literal.PlaintextModulus = 0x101
	params, err := NewParametersFromLiteral(literal)
	require.NoError(t, err)

	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	encoder := NewEncoder(params)
	pt := NewPlaintext(params, 0)
	values := make([]uint64, params.MaxSlots())
	values[0] = 42
	require.NoError(t, encoder.Encode(values, pt))

	ct, err := NewEncryptor(params, sk).EncryptNew(pt)
	require.NoError(t, err)
	nonzeroC1 := false
	for _, row := range ct.Value[1].Coeffs {
		for _, coefficient := range row {
			nonzeroC1 = nonzeroC1 || coefficient != 0
		}
	}
	require.True(t, nonzeroC1, "BGV must retain ordinary RLWE encryption semantics")

	decoded := make([]uint64, len(values))
	require.NoError(t, encoder.Decode(rlwe.NewDecryptor(params, sk).DecryptNew(ct), decoded))
	require.Equal(t, values, decoded)
}
