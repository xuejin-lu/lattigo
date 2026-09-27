package fast

import (
	"math/big"
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

func TestMulIntegerQPrefixRowsExplicitWidthAndValidation(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	in := NewCiphertext(params, 1, 3)
	in.IsNTT = true
	in.Scale = rlwe.NewScale(17)
	for component := range in.Value {
		for row := 0; row < 4; row++ {
			q := params.RingQ().SubRings[row].Modulus
			for i := range in.Value[component].Coeffs[row] {
				in.Value[component].Coeffs[row][i] = (uint64(31+component*19+row*23+i) % q)
			}
		}
	}

	// A three-row operation must leave q3 untouched, including in-place use.
	threeRows := in.CopyNew()
	beforeQ3 := append([]uint64(nil), threeRows.Value[0].Coeffs[3]...)
	wantThree := threeRows.CopyNew()
	require.NoError(t, eval.MulIntegerQPrefixRows(threeRows, big.NewInt(7), 3, threeRows))
	for component := range threeRows.Value {
		for row := 0; row < 3; row++ {
			q := params.RingQ().SubRings[row].Modulus
			for i, value := range wantThree.Value[component].Coeffs[row] {
				want := new(big.Int).Mul(new(big.Int).SetUint64(value), big.NewInt(7))
				want.Mod(want, new(big.Int).SetUint64(q))
				require.Equal(t, want.Uint64(), threeRows.Value[component].Coeffs[row][i])
			}
		}
	}
	require.Equal(t, beforeQ3, threeRows.Value[0].Coeffs[3])
	require.True(t, in.Scale.Equal(threeRows.Scale), "integer multiplication must preserve Scale")

	// Four rows are explicitly supported even when the legacy wrapper would
	// select fewer rows for these parameters.
	fourRows := in.CopyNew()
	wantFour := fourRows.CopyNew()
	require.NoError(t, eval.MulIntegerQPrefixRows(fourRows, big.NewInt(11), 4, fourRows))
	for component := range fourRows.Value {
		for row := 0; row < 4; row++ {
			q := params.RingQ().SubRings[row].Modulus
			for i, value := range wantFour.Value[component].Coeffs[row] {
				want := new(big.Int).Mul(new(big.Int).SetUint64(value), big.NewInt(11))
				want.Mod(want, new(big.Int).SetUint64(q))
				require.Equal(t, want.Uint64(), fourRows.Value[component].Coeffs[row][i])
			}
		}
	}

	beforeInvalid := fourRows.CopyNew()
	require.Error(t, eval.MulIntegerQPrefixRows(fourRows, big.NewInt(3), 5, fourRows))
	for component := range fourRows.Value {
		for row := 0; row < 4; row++ {
			require.Equal(t, beforeInvalid.Value[component].Coeffs[row], fourRows.Value[component].Coeffs[row])
		}
	}
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
