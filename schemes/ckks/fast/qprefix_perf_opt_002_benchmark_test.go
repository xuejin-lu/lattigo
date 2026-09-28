package fast

import (
	"fmt"
	"math/bits"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var qprefixPerfOptUint64Sink uint64
var qprefixPerfOptUint192Sink uint192

func qprefixPerfOptP93Parameters(tb testing.TB) ckks.Parameters {
	tb.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56}, LogDefaultScale: 30,
	})
	require.NoError(tb, err)
	return params
}

type fixedWidthPhaseFixture struct {
	params              ckks.Parameters
	eval                *Evaluator
	input, output       *rlwe.Ciphertext
	rows, targetRows    int
	level, rescaleSteps int
	scratch             fastRescaleScratch
	coeff, staged       []ring.Poly
	restored            []ring.Poly
	magnitudes          [][]uint192
	negative            [][]bool
}

func newFixedWidthPhaseFixture(tb testing.TB, rows int) *fixedWidthPhaseFixture {
	tb.Helper()
	params := qprefixPerfOptP93Parameters(tb)
	level := params.MaxLevel()
	values := rescaleOracleBoundaryValues(tb, params, level, rows)
	input, standardInput := makeRescaleOracleInputs(params, level, rows, values)
	eval := NewEvaluator(params)
	steps := params.LevelsConsumedPerRescaling()
	targetLevel := level - steps
	targetRows := min(rows, qPrefixWidthOrPanic(targetLevel))
	output := NewCiphertext(params, input.Degree(), targetLevel)
	output.IsNTT, output.IsMontgomery = input.IsNTT, input.IsMontgomery
	require.NoError(tb, eval.RescaleQPrefixRows(input, rows, output))
	standardOutput := ckks.NewCiphertext(params, standardInput.Degree(), targetLevel)
	standardOutput.IsNTT, standardOutput.IsMontgomery = standardInput.IsNTT, standardInput.IsMontgomery
	standardEval := ckks.NewEvaluator(params, nil)
	require.NoError(tb, standardEval.Rescale(standardInput, standardOutput))
	fixture := &fixedWidthPhaseFixture{
		params: params, eval: eval, input: input, output: output, rows: rows, targetRows: targetRows,
		level: level, rescaleSteps: steps, scratch: newFastRescaleScratch(params.RingQ()),
		coeff: make([]ring.Poly, len(input.Value)), staged: make([]ring.Poly, len(input.Value)),
		restored: make([]ring.Poly, len(input.Value)), magnitudes: make([][]uint192, len(input.Value)),
		negative: make([][]bool, len(input.Value)),
	}
	for component := range input.Value {
		fixture.coeff[component] = ring.NewPoly(params.N(), rows-1)
		fixture.staged[component] = ring.NewPoly(params.N(), targetRows-1)
		fixture.restored[component] = ring.NewPoly(params.N(), targetRows-1)
		fixture.magnitudes[component] = make([]uint192, params.N())
		fixture.negative[component] = make([]bool, params.N())
	}
	require.NoError(tb, fixture.prefixToCoefficient())
	require.NoError(tb, fixture.fixedWidth())
	fixture.residueStaging()
	fixture.restore()
	for component := range output.Value {
		for row := 0; row < targetRows; row++ {
			require.Equal(tb, output.Value[component].Coeffs[row], fixture.restored[component].Coeffs[row], "Fast fixture component=%d q%d", component, row)
			require.Equal(tb, standardOutput.Value[component].Coeffs[row], output.Value[component].Coeffs[row], "Standard fixture component=%d q%d", component, row)
		}
	}
	return fixture
}

func (f *fixedWidthPhaseFixture) prefixToCoefficient() error {
	for component := range f.input.Value {
		if err := prefixToCoefficientRows(f.params.RingQ(), f.input.Value[component], f.rows, true, false, f.coeff[component]); err != nil {
			return err
		}
	}
	return nil
}

func (f *fixedWidthPhaseFixture) fixedWidth() error {
	for component := range f.input.Value {
		for coefficient := 0; coefficient < f.params.N(); coefficient++ {
			var residues [MaxQPrefixWidth]uint64
			for row := 0; row < f.rows; row++ {
				residues[row] = f.coeff[component].Coeffs[row][coefficient]
			}
			magnitude, negative := centeredQPrefix(reconstructQPrefix(f.rows, residues, &f.scratch), f.scratch.modulus[f.rows-1], f.scratch.half[f.rows-1])
			for step := 0; step < f.rescaleSteps; step++ {
				magnitude = roundedMagnitude192(magnitude, f.params.RingQ().SubRings[f.level-step].Modulus)
				targetWidth, err := QPrefixWidth(f.level - step - 1)
				if err != nil {
					return err
				}
				targetWidth = min(f.rows, targetWidth)
				if cmp192(magnitude, f.scratch.half[targetWidth-1]) > 0 {
					return fmt.Errorf("phase fixture exceeds target capacity at level %d component %d coefficient %d", f.level-step-1, component, coefficient)
				}
			}
			f.magnitudes[component][coefficient] = magnitude
			f.negative[component][coefficient] = negative
		}
	}
	return nil
}

func (f *fixedWidthPhaseFixture) residueStaging() {
	for component := range f.input.Value {
		for coefficient := 0; coefficient < f.params.N(); coefficient++ {
			for row := 0; row < f.targetRows; row++ {
				f.staged[component].Coeffs[row][coefficient] = signedResidue192(f.magnitudes[component][coefficient], f.negative[component][coefficient], f.scratch.q[row])
			}
		}
	}
}

func (f *fixedWidthPhaseFixture) restore() {
	for component := range f.input.Value {
		for row := 0; row < f.targetRows; row++ {
			f.params.RingQ().SubRings[row].NTT(f.staged[component].Coeffs[row], f.restored[component].Coeffs[row])
		}
	}
}

func TestFastRescaleFixedWidthPhaseFixtureMatchesProductionAndStandard(t *testing.T) {
	for _, rows := range []int{2, 4} {
		t.Run(fmt.Sprintf("rows%d", rows), func(t *testing.T) { newFixedWidthPhaseFixture(t, rows) })
	}
}

func BenchmarkFastRescaleWidthCounterfactualLogN13P93(b *testing.B) {
	for _, rows := range []int{2, 4} {
		b.Run(fmt.Sprintf("rows%d", rows), func(b *testing.B) {
			f := newFixedWidthPhaseFixture(b, rows)
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				if err := f.eval.RescaleQPrefixRows(f.input, rows, f.output); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

func BenchmarkFastRescaleFixedWidthPhaseLogN13P93(b *testing.B) {
	for _, rows := range []int{2, 4} {
		b.Run(fmt.Sprintf("rows%d", rows), func(b *testing.B) {
			f := newFixedWidthPhaseFixture(b, rows)
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				if err := f.fixedWidth(); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

func BenchmarkFastQPrefixCRTReconstructionLogN13P93(b *testing.B) {
	params := qprefixPerfOptP93Parameters(b)
	scratch := newFastRescaleScratch(params.RingQ())
	for _, rows := range []int{3, 4} {
		var residues [MaxQPrefixWidth]uint64
		for row := 0; row < rows; row++ {
			residues[row] = (1234567 + uint64(row)*76543) % params.Q()[row]
		}
		b.Run(fmt.Sprintf("rows%d", rows), func(b *testing.B) {
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				qprefixPerfOptUint192Sink = reconstructQPrefix(rows, residues, &scratch)
			}
		})
	}
}

func BenchmarkFastFixedWidthModReductionLogN13P93(b *testing.B) {
	params := qprefixPerfOptP93Parameters(b)
	scratch := newFastRescaleScratch(params.RingQ())
	q := params.Q()
	b.Run("mod128By64", func(b *testing.B) {
		value := scratch.modulus[1]
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			qprefixPerfOptUint64Sink = mod128By64(value.lo, value.mid, q[2])
		}
	})
	b.Run("mod192By64", func(b *testing.B) {
		value := scratch.modulus[2]
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			qprefixPerfOptUint64Sink = mod192By64(value, q[3])
		}
	})
}

func BenchmarkFastRoundedMagnitude192HighLimbDivModLogN13P93(b *testing.B) {
	params := qprefixPerfOptP93Parameters(b)
	q := params.Q()
	modulus, half, _, overflow := q0123Modulus(q[0], q[1], q[2], q[3])
	require.False(b, overflow)
	value := half
	divisor := q[params.MaxLevel()]
	b.Run("source_slash_percent", func(b *testing.B) {
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			qprefixPerfOptUint192Sink = roundedMagnitude192(value, divisor)
		}
	})
	b.Run("exact_bits_Div64_candidate", func(b *testing.B) {
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			qprefixPerfOptUint192Sink = roundedMagnitude192BitsDiv64Candidate(value, divisor)
		}
	})
	_ = modulus
}

func roundedMagnitude192BitsDiv64Candidate(value uint192, divisor uint64) uint192 {
	q2, remainder := bits.Div64(0, value.hi, divisor)
	q1, remainder := bits.Div64(remainder, value.mid, divisor)
	q0, remainder := bits.Div64(remainder, value.lo, divisor)
	quotient := uint192{lo: q0, mid: q1, hi: q2}
	if remainder > divisor/2 {
		var carry uint64
		quotient.lo, carry = bits.Add64(quotient.lo, 1, 0)
		quotient.mid, carry = bits.Add64(quotient.mid, 0, carry)
		quotient.hi, _ = bits.Add64(quotient.hi, 0, carry)
	}
	return quotient
}
