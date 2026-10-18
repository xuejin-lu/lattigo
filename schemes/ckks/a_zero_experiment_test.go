package ckks

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func TestAZeroCKKSMulRelin(t *testing.T) {
	tc := NewTestContext(testInsecurePrec45)
	values0, _, ct0 := tc.NewTestVector(-0.4-0.2i, 0.4+0.2i)
	values1, _, ct1 := tc.NewTestVector(-0.3-0.1i, 0.3+0.1i)

	expected := make([]*bignum.Complex, len(values0))
	mul := bignum.NewComplexMultiplier()
	for i := range values0 {
		expected[i] = bignum.NewComplex()
		mul.Mul(values0[i], values1[i], expected[i])
	}

	require.NoError(t, tc.Evl.MulRelin(ct0, ct1, ct0))
	VerifyTestVectors(tc.Params, tc.Ecd, tc.Dec, expected, ct0, tc.Params.LogDefaultScale(), 0, false, t)
}
