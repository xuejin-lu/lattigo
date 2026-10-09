package fastcore

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func automorphismTestRing(t *testing.T, level int) *ring.Ring {
	t.Helper()
	const n = 16
	generator := ring.NewNTTFriendlyPrimesGenerator(50, 2*n)
	moduli, err := generator.NextAlternatingPrimes(level + 1)
	require.NoError(t, err)
	ringQ, err := ring.NewRing(n, moduli)
	require.NoError(t, err)
	return ringQ.AtLevel(level)
}

func fillAutomorphismTestPoly(ringQ *ring.Ring, poly ring.Poly, seed uint64) {
	for row := range poly.Coeffs {
		q := ringQ.SubRings[row].Modulus
		for i := range poly.Coeffs[row] {
			poly.Coeffs[row][i] = (seed + uint64((row+3)*(i+5))) % q
		}
	}
}

func TestAutomorphismWorkspaceSupportsFullActiveQAndAliases(t *testing.T) {
	const level = 5 // six active rows, beyond the explicit Fast Q-prefix cap
	ringQ := automorphismTestRing(t, level)
	galEl := ringQ.NthRoot() - 1
	rows := level + 1

	for _, isNTT := range []bool{false, true} {
		workspace := NewAutomorphismWorkspace(ringQ.N(), rows)
		in0, in1 := ring.NewPoly(ringQ.N(), level), ring.NewPoly(ringQ.N(), level)
		out0, out1 := ring.NewPoly(ringQ.N(), level), ring.NewPoly(ringQ.N(), level)
		fillAutomorphismTestPoly(ringQ, in0, 17)
		fillAutomorphismTestPoly(ringQ, in1, 41)
		fillAutomorphismTestPoly(ringQ, out0, 73)
		fillAutomorphismTestPoly(ringQ, out1, 89)
		if isNTT {
			for row := 0; row < rows; row++ {
				ringQ.SubRings[row].NTT(in0.Coeffs[row], in0.Coeffs[row])
				ringQ.SubRings[row].NTT(in1.Coeffs[row], in1.Coeffs[row])
			}
		}

		want0, want1 := ring.NewPoly(ringQ.N(), level), ring.NewPoly(ringQ.N(), level)
		if isNTT {
			index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galEl)
			require.NoError(t, err)
			ringQ.AutomorphismNTTWithIndex(in0, index, want0)
			ringQ.AutomorphismNTTWithIndex(in1, index, want1)
		} else {
			ringQ.Automorphism(in0, galEl, want0)
			ringQ.Automorphism(in1, galEl, want1)
		}

		err := workspace.ApplyRows(ringQ, galEl, isNTT, rows,
			PolynomialPair{Input: in0, Output: out0},
			PolynomialPair{Input: in1, Output: out1},
		)
		require.NoError(t, err)
		require.Equal(t, want0.Coeffs, out0.Coeffs)
		require.Equal(t, want1.Coeffs, out1.Coeffs)
		if isNTT {
			require.Equal(t, 1, workspace.CachedGaloisElementCount())
		}

		alias := in0.CopyNew()
		err = workspace.ApplyRows(ringQ, galEl, isNTT, rows, PolynomialPair{Input: *alias, Output: *alias})
		require.NoError(t, err)
		require.Equal(t, want0.Coeffs, alias.Coeffs)
	}
}

func TestAutomorphismWorkspaceValidatesAllPairsBeforeWrites(t *testing.T) {
	const level = 5
	ringQ := automorphismTestRing(t, level)
	rows := level + 1
	workspace := NewAutomorphismWorkspace(ringQ.N(), rows)
	input0, output0 := ring.NewPoly(ringQ.N(), level), ring.NewPoly(ringQ.N(), level)
	input1, output1 := ring.NewPoly(ringQ.N(), level), ring.NewPoly(ringQ.N(), level)
	fillAutomorphismTestPoly(ringQ, input0, 3)
	fillAutomorphismTestPoly(ringQ, input1, 7)
	fillAutomorphismTestPoly(ringQ, output0, 29)
	fillAutomorphismTestPoly(ringQ, output1, 31)
	output0Before := output0.CopyNew()
	output1.Coeffs[level] = nil

	err := workspace.ApplyRows(ringQ, ringQ.NthRoot()-1, false, rows,
		PolynomialPair{Input: input0, Output: output0},
		PolynomialPair{Input: input1, Output: output1},
	)
	require.Error(t, err)
	require.Equal(t, output0Before.Coeffs, output0.Coeffs, "earlier valid pairs must remain unchanged if a later pair is malformed")
}
