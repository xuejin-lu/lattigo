package bootstrapping

import (
	"testing"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func BenchmarkFastPackingLogN13Count2(b *testing.B) {
	benchmarkFastPacking(b, 13, 2)
}

func BenchmarkFastPackingLogN13Count4(b *testing.B) {
	benchmarkFastPacking(b, 13, 4)
}

func BenchmarkFastPackingFactorTwoLogN13Count4(b *testing.B) {
	benchmarkFastPackingFactorTwo(b, 4)
}

func benchmarkFastPacking(b *testing.B, logN, count int) {
	b.Helper()
	generator := ring.NewNTTFriendlyPrimesGenerator(50, uint64(4*(1<<logN)))
	moduli, err := generator.NextAlternatingPrimes(9)
	if err != nil {
		b.Fatal(err)
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: logN, Q: moduli, LogDefaultScale: 30})
	if err != nil {
		b.Fatal(err)
	}
	bootstrapParams := Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: dft.MatrixLiteral{LogSlots: logN - 1},
	}
	eval, err := NewFastEvaluator(bootstrapParams)
	if err != nil {
		b.Fatal(err)
	}
	inputs := make([]rlwe.Ciphertext, count)
	for i := range inputs {
		inputs[i] = *newFastPackingCiphertext(params, 8, logN-3, uint64(601+i*17))
	}
	if _, _, _, err = eval.PackAndSwitchN1ToN2(inputs); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, _, err = eval.PackAndSwitchN1ToN2(inputs); err != nil {
			b.Fatal(err)
		}
	}
}

func benchmarkFastPackingFactorTwo(b *testing.B, count int) {
	b.Helper()
	logN2 := 13
	generator := ring.NewNTTFriendlyPrimesGenerator(50, uint64(4*(1<<logN2)))
	moduli, err := generator.NextAlternatingPrimes(9)
	if err != nil {
		b.Fatal(err)
	}
	n1, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: logN2 - 1, Q: moduli, LogDefaultScale: 30})
	if err != nil {
		b.Fatal(err)
	}
	n2, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: logN2, Q: moduli, LogDefaultScale: 30})
	if err != nil {
		b.Fatal(err)
	}
	bootstrapParams := Parameters{
		ResidualParameters:      n1,
		BootstrappingParameters: n2,
		SlotsToCoeffsParameters: dft.MatrixLiteral{LogSlots: logN2 - 1},
	}
	eval, err := NewFastEvaluator(bootstrapParams)
	if err != nil {
		b.Fatal(err)
	}
	inputs := make([]rlwe.Ciphertext, count)
	for i := range inputs {
		inputs[i] = *newFastPackingCiphertext(n1, 8, logN2-3, uint64(701+i*19))
	}
	if _, _, _, err = eval.PackAndSwitchN1ToN2(inputs); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, _, err = eval.PackAndSwitchN1ToN2(inputs); err != nil {
			b.Fatal(err)
		}
	}
}
