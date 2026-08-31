package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func TestFastCKKSRNSPhase1A(t *testing.T) {
	// Setup a standard RingQ environment with 3 active limbs
	N := 16
	moduli := []uint64{
		0x80000000080001, // q0
		0x2000000a0001,   // q1
		0x2000000e0001,   // q2
	}
	ringQ, err := ring.NewRing(N, moduli)
	require.NoError(t, err)

	backend := NewRNSBackend(ringQ, 2)

	q0 := ringQ.SubRings[0].Modulus
	q1 := ringQ.SubRings[1].Modulus
	Q01 := new(big.Int).Mul(new(big.Int).SetUint64(q0), new(big.Int).SetUint64(q1))
	Q01Half := new(big.Int).Rsh(Q01, 1)

	t.Run("Test 1 - Reconstruction", func(t *testing.T) {
		src := ringQ.NewPoly()

		// Generate zero, positive, negative values across all N coefficients
		expected := make([]*big.Int, N)
		for k := 0; k < N; k++ {
			if k%3 == 0 {
				expected[k] = big.NewInt(0) // Zero
			} else if k%3 == 1 {
				expected[k] = big.NewInt(123456789) // Positive
			} else {
				expected[k] = big.NewInt(-987654321) // Negative
			}

			// Set q0 and q1 residues
			tmp := new(big.Int).Mod(expected[k], new(big.Int).SetUint64(q0))
			if tmp.Sign() < 0 {
				tmp.Add(tmp, new(big.Int).SetUint64(q0))
			}
			src.Coeffs[0][k] = tmp.Uint64()

			tmp = new(big.Int).Mod(expected[k], new(big.Int).SetUint64(q1))
			if tmp.Sign() < 0 {
				tmp.Add(tmp, new(big.Int).SetUint64(q1))
			}
			src.Coeffs[1][k] = tmp.Uint64()
		}

		res, err := backend.ReconstructQ0Q1(&src)
		require.NoError(t, err)

		for k := 0; k < N; k++ {
			require.Equal(t, expected[k].String(), res[k].String())
		}
	})

	t.Run("Test 2 - Only q0/q1 are used", func(t *testing.T) {
		src := ringQ.NewPoly()
		expected := big.NewInt(9999)

		for k := 0; k < N; k++ {
			tmp := new(big.Int).Mod(expected, new(big.Int).SetUint64(q0))
			src.Coeffs[0][k] = tmp.Uint64()
			tmp = new(big.Int).Mod(expected, new(big.Int).SetUint64(q1))
			src.Coeffs[1][k] = tmp.Uint64()

			// Arbitrarily modify q2 to prove reconstruction completely ignores it
			src.Coeffs[2][k] = 123456789 + uint64(k)
		}

		res, err := backend.ReconstructQ0Q1(&src)
		require.NoError(t, err)
		for k := 0; k < N; k++ {
			require.Equal(t, expected.String(), res[k].String())
		}
	})

	t.Run("Test 3 - Centered reconstruction semantics", func(t *testing.T) {
		src := ringQ.NewPoly()

		// We want to test boundary values near Q01/2 and -Q01/2
		posVal := new(big.Int).Sub(Q01Half, big.NewInt(1))
		negVal := new(big.Int).Neg(Q01Half)

		// Set for posVal
		tmp1 := new(big.Int).Mod(posVal, new(big.Int).SetUint64(q0))
		if tmp1.Sign() < 0 {
			tmp1.Add(tmp1, new(big.Int).SetUint64(q0))
		}
		src.Coeffs[0][0] = tmp1.Uint64()

		tmp1 = new(big.Int).Mod(posVal, new(big.Int).SetUint64(q1))
		if tmp1.Sign() < 0 {
			tmp1.Add(tmp1, new(big.Int).SetUint64(q1))
		}
		src.Coeffs[1][0] = tmp1.Uint64()

		// Set for negVal
		tmp2 := new(big.Int).Mod(negVal, new(big.Int).SetUint64(q0))
		if tmp2.Sign() < 0 {
			tmp2.Add(tmp2, new(big.Int).SetUint64(q0))
		}
		src.Coeffs[0][1] = tmp2.Uint64()

		tmp2 = new(big.Int).Mod(negVal, new(big.Int).SetUint64(q1))
		if tmp2.Sign() < 0 {
			tmp2.Add(tmp2, new(big.Int).SetUint64(q1))
		}
		src.Coeffs[1][1] = tmp2.Uint64()

		res, err := backend.ReconstructQ0Q1(&src)
		require.NoError(t, err)

		// Validates that residues correctly map to the unique centered representatives
		require.Equal(t, posVal.String(), res[0].String())
		require.Equal(t, negVal.String(), res[1].String())

		// Test explicit caller-provided boundary checking API on the reconstructed representation
		bound := big.NewInt(1000)
		err = backend.CheckCoefficientBounds(posVal, bound)
		require.Error(t, err) // posVal is near Q01/2 which is >> 1000

		err = backend.CheckCoefficientBounds(big.NewInt(500), bound)
		require.NoError(t, err) // 500 is within 1000
	})

	t.Run("Test 4 - Redistribution validation", func(t *testing.T) {
		x := big.NewInt(-8888)
		coeffs := make([]*big.Int, N)
		for k := 0; k < N; k++ {
			coeffs[k] = new(big.Int).Set(x)
		}

		dst := ringQ.NewPoly()
		err := backend.Redistribute(coeffs, &dst)
		require.NoError(t, err)

		// Validates ALL active limbs correctly inherit residues matching x across all coefficients
		for i := 0; i <= 2; i++ {
			qi := new(big.Int).SetUint64(ringQ.SubRings[i].Modulus)
			tmp := new(big.Int).Mod(x, qi)
			if tmp.Sign() < 0 {
				tmp.Add(tmp, qi)
			}

			for k := 0; k < N; k++ {
				require.Equal(t, tmp.Uint64(), dst.Coeffs[i][k])
			}
		}

		// Error path checks: prevent silent panics
		badCoeffs := make([]*big.Int, N-1) // Length shorter than N
		err = backend.Redistribute(badCoeffs, &dst)
		require.Error(t, err)
		require.Contains(t, err.Error(), "insufficient number of coefficients")

		badDst := ringQ.NewPoly()
		badDst.Coeffs = badDst.Coeffs[:1] // Not enough limbs
		err = backend.Redistribute(coeffs, &badDst)
		require.Error(t, err)
		require.Contains(t, err.Error(), "insufficient RNS limbs")
	})

	t.Run("Test 5 - Full round trip", func(t *testing.T) {
		src := ringQ.NewPoly()
		expectedVals := make([]*big.Int, N)

		for k := 0; k < N; k++ {
			expectedVals[k] = big.NewInt(int64(123456 + k))

			// Fill all original limbs to verify the entire array reconstructs correctly
			// and identically maps back.
			for i := 0; i <= 2; i++ {
				qi := new(big.Int).SetUint64(ringQ.SubRings[i].Modulus)
				tmp := new(big.Int).Mod(expectedVals[k], qi)
				if tmp.Sign() < 0 {
					tmp.Add(tmp, qi)
				}
				src.Coeffs[i][k] = tmp.Uint64()
			}
		}

		dst := ringQ.NewPoly()
		err := backend.ReconstructAndRedistribute(&src, &dst)
		require.NoError(t, err)

		for i := 0; i <= 2; i++ {
			for k := 0; k < N; k++ {
				require.Equal(t, src.Coeffs[i][k], dst.Coeffs[i][k])
			}
		}
	})

	t.Run("Test 6 - Invalid backend configuration", func(t *testing.T) {
		src := ringQ.NewPoly()
		emptyCoeffs := make([]*big.Int, N)

		// 1. NewRNSBackend(nil, ...) checks
		backendNil := NewRNSBackend(nil, 2)
		_, err := backendNil.ReconstructQ0Q1(&src)
		require.Error(t, err)
		require.Contains(t, err.Error(), "RingQ cannot be nil")

		err = backendNil.Redistribute(emptyCoeffs, &src)
		require.Error(t, err)
		require.Contains(t, err.Error(), "RingQ cannot be nil")

		// 2. backend.Level < 0
		backendLow := NewRNSBackend(ringQ, -1)
		_, err = backendLow.ReconstructQ0Q1(&src)
		require.Error(t, err)
		require.Contains(t, err.Error(), "requires at least q0 and q1 (Level >= 1)")

		err = backendLow.Redistribute(emptyCoeffs, &src)
		require.Error(t, err)
		require.Contains(t, err.Error(), "invalid level for redistribution")

		// 3. backend.Level > ringQ.Level()
		backendHigh := NewRNSBackend(ringQ, ringQ.Level()+1)
		_, err = backendHigh.ReconstructQ0Q1(&src)
		require.Error(t, err)
		require.Contains(t, err.Error(), "invalid level: exceeds RingQ.Level()")

		err = backendHigh.Redistribute(emptyCoeffs, &src)
		require.Error(t, err)
		require.Contains(t, err.Error(), "invalid level for redistribution")
	})
}
