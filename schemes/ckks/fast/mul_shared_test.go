package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

type recordingFastMulCore struct {
	delegate fastcore.MulCore
	rows     int
	relin    bool
	calls    int
}

func (recorder *recordingFastMulCore) ApplyRows(ringQ *ring.Ring, op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], opOut *rlwe.Ciphertext, rows int, relin bool) error {
	recorder.rows, recorder.relin = rows, relin
	recorder.calls++
	return recorder.delegate.ApplyRows(ringQ, op0, op1, opOut, rows, relin)
}

func TestExplicitFastMulRelinUsesSharedQPrefixCore(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	recorder := &recordingFastMulCore{delegate: fastcore.NewMulWorkspace(params.RingQ(), MaxQPrefixWidth)}
	eval.mulCore = recorder
	a, b, out := NewCiphertext(params, 1, 3), NewCiphertext(params, 1, 3), NewCiphertext(params, 1, 3)
	a.IsNTT, b.IsNTT, out.IsNTT = true, true, true
	a.Scale, b.Scale = rlwe.NewScale(17), rlwe.NewScale(19)
	for _, ct := range []*rlwe.Ciphertext{a, b} {
		for component := range ct.Value {
			for row := 0; row < 4; row++ {
				q := params.RingQ().SubRings[row].Modulus
				for i := range ct.Value[component].Coeffs[row] {
					ct.Value[component].Coeffs[row][i] = uint64(13+row+component+i) % q
				}
			}
		}
	}
	require.NoError(t, eval.MulRelinElementQPrefixRows(a, b.El(), 4, out))
	require.Equal(t, 1, recorder.calls)
	require.Equal(t, 4, recorder.rows)
	require.True(t, recorder.relin)
	require.Equal(t, 1, out.Degree())
	require.Equal(t, 3, out.Level())
	require.True(t, out.Scale.Equal(a.Scale.Mul(b.Scale)))
}
