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

func BenchmarkFastBootstrapManyLogN13P93Count1(b *testing.B) {
	benchmarkFastBootstrapManyLogN13P93(b, 1)
}

func BenchmarkFastBootstrapManyLogN13P93Count3(b *testing.B) {
	benchmarkFastBootstrapManyLogN13P93(b, 3)
}

func benchmarkFastBootstrapManyLogN13P93(b *testing.B, count int) {
	params, residual := fastLogN13CompressionParameters(b)
	eval, err := NewFastEvaluator(params)
	if err != nil {
		b.Fatal(err)
	}
	inputs := make([]rlwe.Ciphertext, count)
	for i := range inputs {
		values := make([]complex128, 1<<params.LogMaxSlots())
		for j := range values {
			values[j] = complex(float64((j+3*i)%17-8)/256, float64((3*j+i)%13-6)/512)
		}
		inputs[i] = *fastBootstrapEncodedCiphertext(b, residual, params.LogMaxSlots(), values)
	}
	if _, err = eval.BootstrapMany(cloneFastPackingInputs(inputs)); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err = eval.BootstrapMany(cloneFastPackingInputs(inputs)); err != nil {
			b.Fatal(err)
		}
	}
	b.StopTimer()
	b.Logf("profile=LogN13/P93/LogSlots12 batch=%d", count)
}
