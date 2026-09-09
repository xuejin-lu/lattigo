package bootstrapping

import (
	"testing"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func BenchmarkFastEvalModLogN13(b *testing.B) {
	benchmarkFastEvalMod(b, 13)
}

func BenchmarkFastEvalModLogN16(b *testing.B) {
	benchmarkFastEvalMod(b, 16)
}

func benchmarkFastEvalMod(b *testing.B, logN int) {
	logQ := []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60}
	if logN == 16 {
		logQ[0], logQ[1] = 54, 38
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            logQ,
		LogDefaultScale: 45,
	})
	if err != nil {
		b.Fatal(err)
	}
	mod1Params, err := mod1.NewParametersFromLiteral(params, mod1.ParametersLiteral{
		LevelQ:          9,
		LogScale:        60,
		Mod1Type:        mod1.CosDiscrete,
		LogMessageRatio: 2,
		K:               16,
		Mod1Degree:      30,
		DoubleAngle:     3,
	})
	if err != nil {
		b.Fatal(err)
	}
	bootstrapParams := Parameters{BootstrappingParameters: params, Mod1ParametersLiteral: mod1.ParametersLiteral{
		LevelQ:          9,
		LogScale:        60,
		Mod1Type:        mod1.CosDiscrete,
		LogMessageRatio: 2,
		K:               16,
		Mod1Degree:      30,
		DoubleAngle:     3,
	}}
	eval, err := NewFastEvaluator(bootstrapParams)
	if err != nil {
		b.Fatal(err)
	}
	input := newFastEvalModInput(params, mod1Params)
	if _, err = eval.EvalMod(input); err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err = eval.EvalMod(input); err != nil {
			b.Fatal(err)
		}
	}
}
