package bootstrapping

import (
	"strconv"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

func BenchmarkFastBootstrapEndToEnd(b *testing.B) {
	for _, logN := range []int{13, 16} {
		b.Run("LogN="+strconv.Itoa(logN), func(b *testing.B) {
			logSlots := 4
			if logN < logSlots+1 {
				logSlots = logN - 1
			}
			params, residual := fastBootstrapParametersProfileAt(b, logN, logSlots, false, 16, 30, 3)
			params.ResidualParameters = residual
			fastEval, err := NewFastEvaluator(params)
			require.NoError(b, err)
			fastInput := fastBootstrapEncodedCiphertext(b, residual, logSlots, []complex128{0.125 + 0.25i, -0.25 + 0.0625i})
			require.NoError(b, fastEval.ensureFastBootstrapCircuit())

			sk := rlwe.NewKeyGenerator(params.BootstrappingParameters).GenSecretKeyNew()
			keys, _, err := params.GenEvaluationKeys(sk)
			require.NoError(b, err)
			standardEval, err := NewEvaluator(params, keys)
			require.NoError(b, err)
			standardInput := fastBootstrapEncodedCiphertext(b, residual, logSlots, []complex128{0.125 + 0.25i, -0.25 + 0.0625i})

			b.Run("Fast", func(b *testing.B) {
				b.ReportAllocs()
				b.Logf("LogN=%d LogSlots=%d input Level=%d output Level=%d K=%d Mod1Degree=%d DoubleAngle=%d", logN, logSlots, fastInput.Level(), fastEval.OutputLevel(), 16, 30, 3)
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					if _, err = fastEval.Bootstrap(fastInput.CopyNew()); err != nil {
						b.Fatal(err)
					}
				}
			})

			b.Run("Standard", func(b *testing.B) {
				b.ReportAllocs()
				b.Logf("LogN=%d LogSlots=%d input Level=%d output Level=%d K=%d Mod1Degree=%d DoubleAngle=%d", logN, logSlots, standardInput.Level(), standardEval.OutputLevel(), 16, 30, 3)
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					if _, err = standardEval.Bootstrap(standardInput.CopyNew()); err != nil {
						b.Fatal(err)
					}
				}
			})
		})
	}
}
