package polynomial

import (
	"testing"

	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func benchmarkFastPolynomialParameters(b *testing.B, logN int) ckks.Parameters {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            []int{50, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35},
		LogDefaultScale: 30,
	})
	if err != nil {
		b.Fatal(err)
	}
	return params
}

func benchmarkFastPolynomial(b *testing.B, logN int) {
	params := benchmarkFastPolynomialParameters(b, logN)
	eval := NewFastEvaluator(params, nil)
	poly := fastPolynomialTestPoly(30)
	input := fastPolynomialTestCiphertext(params, 29)
	if _, err := eval.Evaluate(input, poly, params.DefaultScale()); err != nil {
		b.Fatal(err)
	}

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := eval.Evaluate(input, poly, params.DefaultScale()); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkFastPolynomialLogN13(b *testing.B) {
	benchmarkFastPolynomial(b, 13)
}

func BenchmarkFastPolynomialLogN16(b *testing.B) {
	benchmarkFastPolynomial(b, 16)
}
