package ckks

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func fastDecryptProjectionParameters(t *testing.T) Parameters {
	t.Helper()
	// Match the BTS LogN16 profile's effective 30-prime Q chain. The decrypt
	// target is the resulting Level-18 ciphertext, but the parameter provider
	// must still describe the complete BTS chain and its later 42-bit limbs.
	logQ := make([]int, 30)
	logQ[0] = 42
	for i := 1; i <= 18; i++ {
		logQ[i] = 32
	}
	for i := 19; i < len(logQ); i++ {
		logQ[i] = 42
	}
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN:            16,
		LogQ:            logQ,
		LogP:            []int{42, 42, 42, 42, 42},
		LogDefaultScale: 32,
	})
	require.NoError(t, err)
	return params
}

func fastDecryptProjectionValues(count int) []float64 {
	values := make([]float64, count)
	for i := range values {
		values[i] = float64(i%31-15) / 64
	}
	return values
}

func newFastDecryptProjectionPlaintext(t *testing.T, params Parameters, level int, values []float64, isNTT bool) *rlwe.Plaintext {
	t.Helper()
	pt := NewPlaintext(params, level)
	pt.IsNTT = isNTT
	require.NoError(t, NewEncoder(params).Encode(values, pt))
	return pt
}

func newFastCompactDecryptProjectionCiphertext(t *testing.T, params Parameters, plaintext *rlwe.Plaintext, logicalLevel, degree int) *rlwe.Ciphertext {
	t.Helper()
	ct := NewCiphertext(params, degree, logicalLevel)
	*ct.MetaData = *plaintext.MetaData
	for component := range ct.Value {
		for row := 0; row <= logicalLevel; row++ {
			if row >= 4 {
				ct.Value[component].Coeffs[row] = nil
			} else if component == 0 {
				copy(ct.Value[component].Coeffs[row], plaintext.Value.Coeffs[row])
			}
		}
	}
	return ct
}

func requireFastDecodedPrefix(t *testing.T, params Parameters, encoder *Encoder, decryptor *rlwe.Decryptor, ct *rlwe.Ciphertext, want []float64, expectedLevel int) *rlwe.Plaintext {
	t.Helper()
	logicalLevel := ct.Level()
	pt := decryptor.DecryptNew(ct)
	require.Equal(t, expectedLevel, pt.Level())
	require.Equal(t, expectedLevel, pt.Value.Level())
	require.Equal(t, logicalLevel, ct.Level(), "projection must not lower the input ciphertext's logical Level")
	require.Equal(t, ct.Scale, pt.Scale)
	require.Equal(t, ct.IsNTT, pt.IsNTT)
	require.Equal(t, ct.IsMontgomery, pt.IsMontgomery)
	require.Equal(t, ct.IsBatched, pt.IsBatched)
	require.Equal(t, ct.LogDimensions, pt.LogDimensions)
	require.Equal(t, ct.IsBitReversed, pt.IsBitReversed)

	rows := min(expectedLevel+1, 4)
	for row := 0; row < rows; row++ {
		require.Equal(t, ct.Value[0].Coeffs[row], pt.Value.Coeffs[row], "coefficient-domain prefix oracle q%d", row)
	}
	for _, component := range ct.Value[1:] {
		for row := 0; row < rows; row++ {
			for _, residue := range component.Coeffs[row] {
				require.Zero(t, residue)
			}
		}
	}

	got := make([]float64, len(want))
	require.NoError(t, encoder.Decode(pt, got))
	for i := range got {
		require.False(t, math.IsNaN(got[i]) || math.IsInf(got[i], 0), "decoded slot %d must be finite", i)
		require.InDelta(t, want[i], got[i], 1e-6, "decoded slot %d", i)
	}
	return pt
}

func TestFastCKKSZeroSecretDecryptNewProjectsCompactLogN16Ciphertext(t *testing.T) {
	params := fastDecryptProjectionParameters(t)
	encoder := NewEncoder(params)
	values := fastDecryptProjectionValues(params.MaxSlots())
	plaintext := newFastDecryptProjectionPlaintext(t, params, 3, values, true)
	sk := rlwe.NewSecretKey(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	ct := newFastCompactDecryptProjectionCiphertext(t, params, plaintext, 18, 2)

	require.Equal(t, 18, ct.Level())
	require.Equal(t, 2, ct.Degree())
	for component := range ct.Value {
		for row := 4; row <= ct.Level(); row++ {
			require.Empty(t, ct.Value[component].Coeffs[row], "dormant q%d must remain absent", row)
		}
	}

	inputBefore := ct.CopyNew()
	decoded := requireFastDecodedPrefix(t, params, encoder, decryptor, ct, values, 3)
	require.True(t, ct.Equal(inputBefore), "DecryptNew must not mutate the ciphertext")

	// The native ciphertext serializer stores zero-length dormant rows, not
	// fabricated residues. A restored snapshot must take the same public
	// DecryptNew -> Decode path and produce the same plaintext and decoded slots.
	snapshot, err := ct.MarshalBinary()
	require.NoError(t, err)
	var restored rlwe.Ciphertext
	require.NoError(t, restored.UnmarshalBinary(snapshot))
	require.Equal(t, 18, restored.Level())
	require.Equal(t, 2, restored.Degree())
	require.True(t, ct.MetaData.Equal(restored.MetaData))
	for component := range restored.Value {
		for row := 0; row < 4; row++ {
			require.Len(t, restored.Value[component].Coeffs[row], params.N())
		}
		for row := 4; row <= restored.Level(); row++ {
			require.Empty(t, restored.Value[component].Coeffs[row], "snapshot must not create dormant q%d residues", row)
		}
	}
	decodedRestored := requireFastDecodedPrefix(t, params, encoder, decryptor, &restored, values, 3)
	require.True(t, decoded.Equal(decodedRestored))
	require.True(t, ct.Equal(inputBefore), "snapshot/restore must not mutate the original ciphertext")
}

func TestFastCKKSZeroSecretDecryptProjectsIntoExistingPlaintext(t *testing.T) {
	params := fastDecryptProjectionParameters(t)
	encoder := NewEncoder(params)
	values := fastDecryptProjectionValues(256)
	plaintext := newFastDecryptProjectionPlaintext(t, params, 3, values, false)
	sk := rlwe.NewSecretKey(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	ct := newFastCompactDecryptProjectionCiphertext(t, params, plaintext, 18, 2)

	out := NewPlaintext(params, 18)
	out.Value.Coeffs[0][0] = 123
	outBefore := out.CopyNew()
	ctBefore := ct.CopyNew()
	decryptor.Decrypt(ct, out)
	require.Equal(t, 3, out.Level())
	require.Equal(t, 3, out.Value.Level())
	require.True(t, out.MetaData.Equal(ct.MetaData))
	require.True(t, ct.Equal(ctBefore))
	got := make([]float64, len(values))
	require.NoError(t, encoder.Decode(out, got))
	for i := range got {
		require.InDelta(t, values[i], got[i], 1e-6)
	}
	require.NotEqual(t, outBefore.Value.Coeffs[0][0], out.Value.Coeffs[0][0])
}

func TestFastCKKSZeroSecretDecryptLowLevelBoundariesAndCompleteQPath(t *testing.T) {
	params := fastDecryptProjectionParameters(t)
	encoder := NewEncoder(params)
	sk := rlwe.NewSecretKey(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	values := fastDecryptProjectionValues(128)

	for level := 0; level <= 3; level++ {
		t.Run(fmt.Sprintf("level_%d", level), func(t *testing.T) {
			plaintext := newFastDecryptProjectionPlaintext(t, params, level, values, level%2 == 0)
			ct := NewCiphertext(params, 1, level)
			*ct.MetaData = *plaintext.MetaData
			copy(ct.Value[0].Coeffs[0], plaintext.Value.Coeffs[0])
			for row := 1; row <= level; row++ {
				copy(ct.Value[0].Coeffs[row], plaintext.Value.Coeffs[row])
			}
			got := decryptor.DecryptNew(ct)
			require.Equal(t, level, got.Level())
			require.Equal(t, level, got.Value.Level())
			decoded := make([]float64, len(values))
			require.NoError(t, encoder.Decode(got, decoded))
			for i := range decoded {
				require.InDelta(t, values[i], decoded[i], 1e-6)
			}
		})
	}

	// A complete logical-Q ciphertext remains on the generic decrypt path and
	// retains its original high plaintext Level even under Fast CKKS params.
	plaintext := newFastDecryptProjectionPlaintext(t, params, 3, values, true)
	full := NewCiphertext(params, 1, 18)
	*full.MetaData = *plaintext.MetaData
	for row := 0; row <= 18; row++ {
		if row < 4 {
			copy(full.Value[0].Coeffs[row], plaintext.Value.Coeffs[row])
		}
	}
	fullPlaintext := decryptor.DecryptNew(full)
	require.Equal(t, 18, fullPlaintext.Level())
	require.Equal(t, 18, fullPlaintext.Value.Level())

	// A core RLWE ParameterProvider has no Fast capability marker and therefore
	// retains the ordinary full-Q decrypt behavior too.
	genericDecryptor := rlwe.NewDecryptor(params.GetRLWEParameters(), sk)
	genericPlaintext := genericDecryptor.DecryptNew(full)
	require.Equal(t, 18, genericPlaintext.Level())
}

func TestFastCKKSZeroSecretDecryptRejectsMalformedCompactInputsTransactionally(t *testing.T) {
	params := fastDecryptProjectionParameters(t)
	values := fastDecryptProjectionValues(64)
	plaintext := newFastDecryptProjectionPlaintext(t, params, 3, values, true)
	sk := rlwe.NewSecretKey(params)
	decryptor := rlwe.NewDecryptor(params, sk)

	tests := []struct {
		name   string
		mutate func(*rlwe.Ciphertext)
	}{
		{name: "nonzero_c1", mutate: func(ct *rlwe.Ciphertext) { ct.Value[1].Coeffs[0][0] = 1 }},
		{name: "nonzero_c2", mutate: func(ct *rlwe.Ciphertext) { ct.Value[2].Coeffs[3][0] = 1 }},
		{name: "missing_authoritative_q2", mutate: func(ct *rlwe.Ciphertext) { ct.Value[0].Coeffs[2] = nil }},
		{name: "noncanonical_c0", mutate: func(ct *rlwe.Ciphertext) { ct.Value[0].Coeffs[1][0] = params.RingQ().SubRings[1].Modulus }},
		{name: "zero_scale", mutate: func(ct *rlwe.Ciphertext) { ct.Scale = rlwe.NewScale(0) }},
		{name: "invalid_slot_dimensions", mutate: func(ct *rlwe.Ciphertext) { ct.LogDimensions.Cols = params.LogN() }},
		{name: "unsupported_montgomery_domain", mutate: func(ct *rlwe.Ciphertext) { ct.IsMontgomery = true }},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			ct := newFastCompactDecryptProjectionCiphertext(t, params, plaintext, 18, 2)
			test.mutate(ct)
			invalidBefore := ct.CopyNew()

			out := NewPlaintext(params, 18)
			out.Value.Coeffs[0][0] = 12345
			outBefore := out.CopyNew()
			require.Panics(t, func() { decryptor.Decrypt(ct, out) })
			require.True(t, ct.Equal(invalidBefore), "failed decrypt must not modify input")
			require.True(t, out.Equal(outBefore), "failed decrypt must not modify existing output")
		})
	}

	ct := newFastCompactDecryptProjectionCiphertext(t, params, plaintext, 18, 2)
	ct.MetaData = nil
	require.Panics(t, func() { decryptor.DecryptNew(ct) }, "nil metadata must be rejected explicitly")
}

func TestFastCKKSZeroSecretDecryptWithKeyRetainsProjectionCapability(t *testing.T) {
	params := fastDecryptProjectionParameters(t)
	values := fastDecryptProjectionValues(32)
	plaintext := newFastDecryptProjectionPlaintext(t, params, 3, values, true)
	sk := rlwe.NewSecretKey(params)
	// This branch must not read the secret key at all: the compact c1 prefix is
	// proven zero, so even a key with only q0 backing is sufficient for the test.
	sk.Value.Q.Coeffs = sk.Value.Q.Coeffs[:1]
	decryptor := rlwe.NewDecryptor(params, sk).WithKey(sk)
	ct := newFastCompactDecryptProjectionCiphertext(t, params, plaintext, 18, 1)
	decoded := decryptor.DecryptNew(ct)
	require.Equal(t, 3, decoded.Level())
}

func TestFastCKKSZeroSecretDecryptPreservesCompleteQRingTypes(t *testing.T) {
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN:            8,
		LogQ:            []int{40, 35, 35, 35, 35},
		RingType:        ring.ConjugateInvariant,
		LogDefaultScale: 32,
	})
	require.NoError(t, err)
	sk := rlwe.NewSecretKey(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	// Complete-Q ciphertexts retain the generic decrypt path for every CKKS
	// ring type; the experimental prefix projection is not involved.
	fullQ := NewCiphertext(params, 1, 4)
	fullQ.Scale = rlwe.NewScale(1)
	fullQ.LogDimensions = ring.Dimensions{Cols: params.LogN()}
	decoded := decryptor.DecryptNew(fullQ)
	require.Equal(t, 4, decoded.Level())

	// The compact projection itself is intentionally restricted to the
	// Standard ring, since the CKKS decoder contract here was validated only
	// for that representation.
	ct := NewCiphertext(params, 1, 4)
	ct.Scale = rlwe.NewScale(1)
	ct.LogDimensions = ring.Dimensions{Cols: params.LogN()}
	for row := 1; row <= 4; row++ {
		ct.Value[0].Coeffs[row] = nil
		ct.Value[1].Coeffs[row] = nil
	}
	require.Panics(t, func() { decryptor.DecryptNew(ct) })
}
