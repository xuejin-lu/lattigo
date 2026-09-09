package bootstrapping

import (
	"strconv"
	"testing"

	"github.com/stretchr/testify/require"
)

func BenchmarkFastBootstrapEndToEnd(b *testing.B) {
	for _, logN := range []int{13, 16} {
		b.Run("LogN="+strconv.Itoa(logN), func(b *testing.B) {
			logSlots := 4
			if logN < logSlots+1 {
				logSlots = logN - 1
			}
			params, residual := fastBootstrapParametersAt(b, logN, logSlots, false)
			params.ResidualParameters = residual
			eval, err := NewFastEvaluator(params)
			require.NoError(b, err)
			ct := fastBootstrapEncodedCiphertext(b, residual, logSlots, []complex128{0.125 + 0.25i, -0.25 + 0.0625i})
			require.NoError(b, eval.ensureFastBootstrapCircuit())
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				if _, err = eval.Bootstrap(ct); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
