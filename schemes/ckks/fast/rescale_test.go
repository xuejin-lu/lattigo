package fast

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func rescaleTestParameters(t *testing.T) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 50, 50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	return params
}

func TestFixedWidthCRTMatchesBigInt(t *testing.T) {
	params := rescaleTestParameters(t)
	q0 := params.Q()[0]
	q1 := params.Q()[1]
	inv, ok := inverseMod(q0%q1, q1)
	require.True(t, ok)
	Q := new(big.Int).Mul(new(big.Int).SetUint64(q0), new(big.Int).SetUint64(q1))
	values := []*big.Int{
		big.NewInt(0), big.NewInt(1), big.NewInt(-1), big.NewInt(7),
		new(big.Int).Rsh(new(big.Int).Set(Q), 1),
		new(big.Int).Sub(new(big.Int).Rsh(new(big.Int).Set(Q), 1), big.NewInt(1)),
		new(big.Int).Neg(new(big.Int).Rsh(new(big.Int).Set(Q), 1)),
		new(big.Int).Add(new(big.Int).Neg(new(big.Int).Rsh(new(big.Int).Set(Q), 1)), big.NewInt(1)),
	}
	for i := 0; i < 128; i++ {
		values = append(values, new(big.Int).Lsh(big.NewInt(int64(i+1)), uint(i%17)))
		values = append(values, new(big.Int).Neg(new(big.Int).Lsh(big.NewInt(int64(i+1)), uint(i%17))))
	}
	for _, want := range values {
		r0 := new(big.Int).Mod(new(big.Int).Set(want), new(big.Int).SetUint64(q0)).Uint64()
		r1 := new(big.Int).Mod(new(big.Int).Set(want), new(big.Int).SetUint64(q1)).Uint64()
		lo, hi := crtQ01(r0, r1, q0, q1, inv)
		got := new(big.Int).SetUint64(hi)
		got.Lsh(got, 64)
		got.Add(got, new(big.Int).SetUint64(lo))
		if got.Cmp(new(big.Int).Rsh(new(big.Int).Set(Q), 1)) > 0 {
			got.Sub(got, Q)
		}
		require.Equal(t, want, got)
	}
}

func TestFixedWidthRoundedDivisionMatchesRule(t *testing.T) {
	divisors := []uint64{paramsModulus(50), paramsModulus(49)}
	q01 := new(big.Int).Mul(new(big.Int).SetUint64(paramsModulus(55)), new(big.Int).SetUint64(paramsModulus(39)))
	q01Lo := q01.Uint64()
	q01Hi := new(big.Int).Rsh(new(big.Int).Set(q01), 64).Uint64()
	half := new(big.Int).Rsh(new(big.Int).Set(q01), 1)
	halfLo := half.Uint64()
	halfHi := new(big.Int).Rsh(new(big.Int).Set(half), 64).Uint64()
	for _, d := range divisors {
		for _, k := range []uint64{0, 1, 3, 17} {
			for _, delta := range []int64{-1, 0, 1} {
				m := k * d
				if delta < 0 {
					m -= uint64(-delta)
				} else {
					m += uint64(delta)
				}
				positive := new(big.Int).SetUint64(m)
				lo, hi := positive.Uint64(), uint64(0)
				got, negative := roundedMagnitude128(lo, hi, q01Lo, q01Hi, halfLo, halfHi, d)
				require.False(t, negative)
				want := m / d
				if m%d > d/2 {
					want++
				}
				require.Equal(t, want, got)

				negativeValue := new(big.Int).Sub(q01, positive)
				got, negative = roundedMagnitude128(negativeValue.Uint64(), new(big.Int).Rsh(new(big.Int).Set(negativeValue), 64).Uint64(), q01Lo, q01Hi, halfLo, halfHi, d)
				require.True(t, negative)
				require.Equal(t, want, got)
			}
		}
	}
}

func TestFastRescaleMatchesStandardAndReachesLevelZero(t *testing.T) {
	params := rescaleTestParameters(t)
	standard := ckks.NewEvaluator(params, nil)
	fastEval := NewEvaluator(params)
	values := []int64{0, 1, -1, 12345, -67890, 1 << 20, -(1 << 20)}
	in := makeFastRescaleCiphertext(params, 1, values, 9)
	standardOut := ckks.NewCiphertext(params, 1, 0)
	fastOut := ckks.NewCiphertext(params, 1, 0)
	require.NoError(t, standard.Rescale(in, standardOut))
	require.NoError(t, fastEval.Rescale(in, fastOut))
	require.Equal(t, standardOut.Scale, fastOut.Scale)
	require.Equal(t, standardOut.Level(), fastOut.Level())
	require.Equal(t, standardOut.Value[0].Coeffs[0], fastOut.Value[0].Coeffs[0])
	require.Equal(t, standardOut.Value[1].Coeffs[0], fastOut.Value[1].Coeffs[0])
	require.True(t, fastOut.IsNTT)
	require.False(t, fastOut.IsMontgomery)

	levelOne := makeFastRescaleCiphertext(params, 1, values, 11)
	levelOne.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 120))
	out := ckks.NewCiphertext(params, 1, 1)
	require.NoError(t, fastEval.RescaleTo(levelOne, rlwe.NewScale(1), out))
	require.Equal(t, 0, out.Level())
}

func TestFastRescaleMatchesStandardAtHigherLevels(t *testing.T) {
	params := rescaleTestParameters(t)
	standard := ckks.NewEvaluator(params, nil)
	fastEval := NewEvaluator(params)
	for _, level := range []int{2, 3} {
		d := params.RingQ().SubRings[level].Modulus
		positive := []uint64{d - 1, d + 1, 2*d + 1, 3*d + (d-1)/2, 3*d + (d+1)/2}
		values := make([]int64, 0, 2*len(positive))
		for _, value := range positive {
			values = append(values, int64(value), -int64(value))
		}
		in := makeFastRescaleCiphertext(params, level, values, 0)
		standardOut := ckks.NewCiphertext(params, 1, level-1)
		fastOut := ckks.NewCiphertext(params, 1, level-1)
		require.NoError(t, standard.Rescale(in, standardOut))
		require.NoError(t, fastEval.Rescale(in, fastOut))
		require.Equal(t, standardOut.Scale, fastOut.Scale)
		for component := range fastOut.Value {
			require.NotEqual(t, make([]uint64, len(fastOut.Value[component].Coeffs[0])), fastOut.Value[component].Coeffs[0])
			for limb := 0; limb <= level-1 && limb < 2; limb++ {
				require.Equal(t, standardOut.Value[component].Coeffs[limb], fastOut.Value[component].Coeffs[limb])
			}
		}
	}
}

func TestFastRescaleToSequentialLevelsAndInPlace(t *testing.T) {
	params := rescaleTestParameters(t)
	standard := ckks.NewEvaluator(params, nil)
	fastEval := NewEvaluator(params)
	in := makeFastRescaleCiphertext(params, 3, []int64{1, -2, 12345, -67890}, 0)
	in.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 160))
	standardOut := ckks.NewCiphertext(params, 1, 3)
	fastOut := ckks.NewCiphertext(params, 1, 3)
	minScale := rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 20))
	require.NoError(t, standard.RescaleTo(in, minScale, standardOut))
	require.NoError(t, fastEval.RescaleTo(in, minScale, fastOut))
	require.Equal(t, standardOut.Level(), fastOut.Level())
	require.Equal(t, standardOut.Scale, fastOut.Scale)
	for component := range fastOut.Value {
		for limb := 0; limb <= fastOut.Level() && limb < 2; limb++ {
			require.Equal(t, standardOut.Value[component].Coeffs[limb], fastOut.Value[component].Coeffs[limb])
		}
	}

	inPlace := in.CopyNew()
	outOfPlace := ckks.NewCiphertext(params, 1, 3)
	require.NoError(t, fastEval.RescaleTo(inPlace, minScale, inPlace))
	require.NoError(t, fastEval.RescaleTo(in, minScale, outOfPlace))
	require.Equal(t, outOfPlace.Level(), inPlace.Level())
	for component := range inPlace.Value {
		for limb := 0; limb <= inPlace.Level() && limb < 2; limb++ {
			require.Equal(t, outOfPlace.Value[component].Coeffs[limb], inPlace.Value[component].Coeffs[limb])
		}
	}
}

func TestFastRescaleSupportsBothNTTRepresentations(t *testing.T) {
	params := rescaleTestParameters(t)
	eval := NewEvaluator(params)
	nonMontgomery := makeFastRescaleCiphertext(params, 2, []int64{1, -2, 12345, -67890}, 41)
	montgomery := makeFastRescaleCiphertext(params, 2, []int64{1, -2, 12345, -67890}, 41, true)
	nonOut := ckks.NewCiphertext(params, 1, 1)
	montOut := ckks.NewCiphertext(params, 1, 1)
	require.NoError(t, eval.Rescale(nonMontgomery, nonOut))
	require.NoError(t, eval.Rescale(montgomery, montOut))
	require.False(t, nonOut.IsMontgomery)
	require.True(t, montOut.IsMontgomery)
	r := params.RingQ().AtLevel(1)
	for component := range nonOut.Value {
		for limb := 0; limb < 2; limb++ {
			converted := append([]uint64(nil), montOut.Value[component].Coeffs[limb]...)
			r.SubRings[limb].IMForm(converted, converted)
			require.Equal(t, nonOut.Value[component].Coeffs[limb], converted)
		}
	}
}

func TestFastRescaleIgnoresDormantResidues(t *testing.T) {
	params := rescaleTestParameters(t)
	fastEval := NewEvaluator(params)
	a := makeFastRescaleCiphertext(params, 2, []int64{1, -2, 3, -4}, 17)
	b := a.CopyNew()
	for component := range b.Value {
		for limb := 2; limb < len(b.Value[component].Coeffs); limb++ {
			for k := range b.Value[component].Coeffs[limb] {
				b.Value[component].Coeffs[limb][k] = uint64(1000 + limb + k)
			}
		}
	}
	outA := ckks.NewCiphertext(params, 1, 1)
	outB := ckks.NewCiphertext(params, 1, 1)
	require.NoError(t, fastEval.Rescale(a, outA))
	require.NoError(t, fastEval.Rescale(b, outB))
	for component := range outA.Value {
		require.Equal(t, outA.Value[component].Coeffs[0], outB.Value[component].Coeffs[0])
		require.Equal(t, outA.Value[component].Coeffs[1], outB.Value[component].Coeffs[1])
	}
}

func makeFastRescaleCiphertext(params ckks.Parameters, level int, values []int64, poison uint64, montgomery ...bool) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, level)
	ct.IsNTT = true
	ct.IsMontgomery = len(montgomery) > 0 && montgomery[0]
	ct.Scale = rlwe.NewScale(1 << 40)
	r := params.RingQ().AtLevel(level)
	for component := range ct.Value {
		for limb := 0; limb <= level; limb++ {
			q := r.SubRings[limb].Modulus
			for k := range ct.Value[component].Coeffs[limb] {
				v := values[(k+component)%len(values)]
				if v < 0 {
					ct.Value[component].Coeffs[limb][k] = q - uint64(-v)%q
				} else {
					ct.Value[component].Coeffs[limb][k] = uint64(v) % q
				}
				if limb >= 2 && poison != 0 {
					ct.Value[component].Coeffs[limb][k] += poison + uint64(limb+k)
					ct.Value[component].Coeffs[limb][k] %= q
				}
			}
		}
		for limb := 0; limb <= level; limb++ {
			ntt := make([]uint64, len(ct.Value[component].Coeffs[limb]))
			r.SubRings[limb].NTT(ct.Value[component].Coeffs[limb], ntt)
			if ct.IsMontgomery {
				r.SubRings[limb].MForm(ntt, ntt)
			}
			copy(ct.Value[component].Coeffs[limb], ntt)
		}
	}
	return ct
}

func paramsModulus(bits int) uint64 {
	return (uint64(1) << uint(bits)) - 1
}

func BenchmarkFastRescale(b *testing.B) {
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: 4, LogQ: []int{55, 39, 50}, LogDefaultScale: 30})
	if err != nil {
		b.Fatal(err)
	}
	eval := NewEvaluator(params)
	ct := makeFastRescaleCiphertext(params, 1, []int64{1, -2, 3, -4}, 7)
	out := ckks.NewCiphertext(params, 1, 0)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if err := eval.Rescale(ct, out); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkFastRescaleLogN13(b *testing.B) {
	benchmarkRescaleLogN(b, 13, true)
}

func BenchmarkStandardRescaleLogN13(b *testing.B) {
	benchmarkRescaleLogN(b, 13, false)
}

func BenchmarkFastRescaleLogN16(b *testing.B) {
	benchmarkRescaleLogN(b, 16, true)
}

func BenchmarkStandardRescaleLogN16(b *testing.B) {
	benchmarkRescaleLogN(b, 16, false)
}

func benchmarkRescaleLogN(b *testing.B, logN int, fastPath bool) {
	logQ := []int{55, 39, 39, 45, 60, 60, 60, 60, 60, 60, 60, 60, 56, 56, 56, 56}
	if logN == 16 {
		// The N=65536 NTT-prime generator can select q0/q1 just above the
		// documented profile; use the nearest profile within the fixed-width
		// production contract for this optional benchmark.
		logQ[0], logQ[1] = 54, 38
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{LogN: logN, LogQ: logQ, LogDefaultScale: 30})
	if err != nil {
		b.Fatal(err)
	}
	level := params.MaxLevel()
	in := makeFastRescaleCiphertext(params, level, []int64{1, -2, 12345, -67890}, 0)
	out := ckks.NewCiphertext(params, 1, level-1)
	if fastPath {
		eval := NewEvaluator(params)
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(in, out); err != nil {
				b.Fatal(err)
			}
		}
	} else {
		eval := ckks.NewEvaluator(params, nil)
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(in, out); err != nil {
				b.Fatal(err)
			}
		}
	}
}
