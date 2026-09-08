package fast

import (
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func testFastCKKSParameters(t *testing.T) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{50, 50, 50, 50},
		LogDefaultScale: 40,
	})
	require.NoError(t, err)
	return params
}

func fillFastCKKSCiphertext(ct *rlwe.Ciphertext, params ckks.Parameters, values [][]int64, offset uint64) {
	r := params.RingQ().AtLevel(ct.Level())
	for d := range ct.Value {
		for i, s := range r.SubRings[:ct.Level()+1] {
			for j := range ct.Value[d].Coeffs[i] {
				if d < len(values) && j < len(values[d]) {
					v := values[d][j]
					if v < 0 {
						ct.Value[d].Coeffs[i][j] = s.Modulus - uint64(-v)
					} else {
						ct.Value[d].Coeffs[i][j] = uint64(v)
					}
				} else {
					ct.Value[d].Coeffs[i][j] = (offset + uint64(13*d+7*i+j)) % s.Modulus
				}
			}
		}
	}
}

func fillFastCKKSCiphertextBenchmark(ct *rlwe.Ciphertext, params ckks.Parameters, offset uint64) {
	for d := range ct.Value {
		for i, s := range params.RingQ().SubRings[:ct.Level()+1] {
			for j := range ct.Value[d].Coeffs[i] {
				ct.Value[d].Coeffs[i][j] = (offset + uint64(j+3*i+17*d)) % s.Modulus
			}
		}
	}
}

func benchmarkFastCKKSParameters(b *testing.B) ckks.Parameters {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 10, LogQ: []int{50, 50, 50, 50, 50, 50, 50, 50, 50}, LogDefaultScale: 40})
	if err != nil {
		b.Fatal(err)
	}
	return params
}
