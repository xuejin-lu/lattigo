package fast

import (
	"math/bits"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// BenchmarkFastRescaleQPrefixRows4LogN13P93 measures the four-row, high-level
// q0=55 Rescale surface used by the generated-power path. Input construction
// and the evaluator warmup are outside the timed region.
func BenchmarkFastRescaleQPrefixRows4LogN13P93(b *testing.B) {
	logQ := []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 13, LogQ: logQ, LogDefaultScale: 30})
	require.NoError(b, err)
	require.Equal(b, 55, bits.Len64(params.Q()[0]), "benchmark profile must use q0=55")
	level := params.MaxLevel()
	values := rescaleOracleBoundaryValues(b, params, level, 4)
	input, _ := makeRescaleOracleInputs(params, level, 4, values)
	eval := NewEvaluator(params)
	output := NewCiphertext(params, input.Degree(), level-1)
	if err = eval.RescaleQPrefixRows(input, 4, output); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if err = eval.RescaleQPrefixRows(input, 4, output); err != nil {
			b.Fatal(err)
		}
	}
}
