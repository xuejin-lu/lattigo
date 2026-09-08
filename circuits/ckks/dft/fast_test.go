package dft

import (
	"fmt"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func fastDFTTestParameters(t testing.TB) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{55, 39, 50, 50, 50},
		LogP:            []int{50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	return params
}

func fastDFTMatrix(t testing.TB, params ckks.Parameters, typ Type, format Format) Matrix {
	t.Helper()
	literal := MatrixLiteral{
		Type:         typ,
		LogSlots:     2,
		LevelQ:       3,
		LevelP:       0,
		Levels:       []int{1, 1},
		Format:       format,
		LogBSGSRatio: 1,
	}
	matrix, err := NewMatrixFromLiteral(params, literal, ckks.NewEncoder(params))
	require.NoError(t, err)
	return matrix
}

func fastDFTInput(params ckks.Parameters, level int) *rlwe.Ciphertext {
	ct := ckks.NewCiphertext(params, 1, level)
	ct.IsNTT, ct.IsMontgomery, ct.IsBatched = true, false, true
	ct.Scale = params.DefaultScale()
	ct.LogDimensions = ring.Dimensions{Cols: 2}
	for limb, subring := range params.RingQ().SubRings[:level+1] {
		for i := range ct.Value[0].Coeffs[limb] {
			ct.Value[0].Coeffs[limb][i] = new(big.Int).Mod(big.NewInt(int64((i%7)-3)), new(big.Int).SetUint64(subring.Modulus)).Uint64()
		}
	}
	ct.Value[1].Zero()
	params.RingQ().AtLevel(level).NTT(ct.Value[0], ct.Value[0])
	params.RingQ().AtLevel(level).NTT(ct.Value[1], ct.Value[1])
	return ct
}

func fastDFTMontgomeryCopy(params ckks.Parameters, src *rlwe.Ciphertext) *rlwe.Ciphertext {
	ct := src.CopyNew()
	for d := range ct.Value {
		for limb := 0; limb <= ct.Level(); limb++ {
			params.RingQ().SubRings[limb].MForm(ct.Value[d].Coeffs[limb], ct.Value[d].Coeffs[limb])
		}
	}
	ct.IsMontgomery = true
	return ct
}

func fastDFTNonMontgomeryCopy(params ckks.Parameters, src *rlwe.Ciphertext) *rlwe.Ciphertext {
	ct := src.CopyNew()
	if ct.IsMontgomery {
		for d := range ct.Value {
			for limb := 0; limb <= ct.Level(); limb++ {
				params.RingQ().SubRings[limb].IMForm(ct.Value[d].Coeffs[limb], ct.Value[d].Coeffs[limb])
			}
		}
		ct.IsMontgomery = false
	}
	return ct
}

func requireFastDFTMatchesStandard(t *testing.T, params ckks.Parameters, fast *rlwe.Ciphertext, standard *rlwe.Ciphertext) {
	t.Helper()
	got := fastDFTNonMontgomeryCopy(params, fast)
	require.Equal(t, standard.Level(), got.Level())
	require.Equal(t, standard.Scale, got.Scale)
	require.Equal(t, standard.LogDimensions, got.LogDimensions)
	for d := range standard.Value {
		require.Equal(t, standard.Value[d].Coeffs[:2], got.Value[d].Coeffs[:2], "component %d", d)
	}
}

func standardDFTEvaluator(t testing.TB, params ckks.Parameters, matrices ...Matrix) *Evaluator {
	t.Helper()
	kgen := rlwe.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	var galEls []uint64
	for _, matrix := range matrices {
		galEls = append(galEls, matrix.GaloisElements(params)...)
	}
	evk := rlwe.NewMemEvaluationKeySet(nil, kgen.GenGaloisKeysNew(galEls, sk)...)
	return NewEvaluator(params, ckks.NewEvaluator(params, evk))
}

func TestFastDFTStandardCoeffsToSlotsAndSlotsToCoeffs(t *testing.T) {
	params := fastDFTTestParameters(t)
	cts := fastDFTMatrix(t, params, HomomorphicEncode, Standard)
	stc := fastDFTMatrix(t, params, HomomorphicDecode, Standard)
	fastEval := NewFastEvaluator(params)
	standardEval := standardDFTEvaluator(t, params, cts, stc)
	standardInput := fastDFTInput(params, cts.LevelQ)
	fastInput := fastDFTMontgomeryCopy(params, standardInput)

	fastReal, fastImag, err := fastEval.CoeffsToSlotsNew(fastInput, cts)
	require.NoError(t, err)
	require.Nil(t, fastImag)
	standardReal, standardImag, err := standardEval.CoeffsToSlotsNew(standardInput, cts)
	require.NoError(t, err)
	require.Nil(t, standardImag)
	requireFastDFTMatchesStandard(t, params, fastReal, standardReal)
	require.Equal(t, cts.LevelQ-len(cts.Levels)*params.LevelsConsumedPerRescaling(), fastReal.Level())

	fastDecoded, err := fastEval.SlotsToCoeffsNew(fastInput, nil, stc)
	require.NoError(t, err)
	standardDecoded, err := standardEval.SlotsToCoeffsNew(standardInput, nil, stc)
	require.NoError(t, err)
	requireFastDFTMatchesStandard(t, params, fastDecoded, standardDecoded)
}

func TestFastDFTSplitRepackAndPoison(t *testing.T) {
	params := fastDFTTestParameters(t)
	matrix := fastDFTMatrix(t, params, HomomorphicEncode, RepackImagAsReal)
	fastEval := NewFastEvaluator(params)
	in := fastDFTMontgomeryCopy(params, fastDFTInput(params, matrix.LevelQ))
	poisoned := in.CopyNew()
	for d := range poisoned.Value {
		for limb := 2; limb <= poisoned.Level(); limb++ {
			for i := range poisoned.Value[d].Coeffs[limb] {
				poisoned.Value[d].Coeffs[limb][i] = ^uint64(0) - uint64(i+13*limb+d)
			}
		}
	}
	fastReal, fastImag, err := fastEval.CoeffsToSlotsNew(in, matrix)
	require.NoError(t, err)
	require.Nil(t, fastImag)
	poisonedReal, poisonedImag, err := fastEval.CoeffsToSlotsNew(poisoned, matrix)
	require.NoError(t, err)
	require.Nil(t, poisonedImag)
	require.Equal(t, fastReal.Value[0].Coeffs[:2], poisonedReal.Value[0].Coeffs[:2])
	require.Equal(t, fastReal.Value[1].Coeffs[:2], poisonedReal.Value[1].Coeffs[:2])
	require.Equal(t, fastReal.Level(), poisonedReal.Level())
	require.Equal(t, fastReal.Scale, poisonedReal.Scale)
}

func BenchmarkFastDFTSequence(b *testing.B) {
	params := fastDFTTestParameters(b)
	matrices := fastDFTMatrix(b, params, HomomorphicEncode, Standard)
	standardEval := standardDFTEvaluator(b, params, matrices)
	fastEval := NewFastEvaluator(params)
	standardInput := fastDFTInput(params, matrices.LevelQ)
	fastInput := fastDFTMontgomeryCopy(params, standardInput)
	b.Run("FastQ01", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, _, err := fastEval.CoeffsToSlotsNew(fastInput, matrices); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("StandardFullRNS", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, _, err := standardEval.CoeffsToSlotsNew(standardInput, matrices); err != nil {
				b.Fatal(err)
			}
		}
	})
}

func fastDFTBenchmarkParameters(t testing.TB, logN int) ckks.Parameters {
	t.Helper()
	logQ := []int{55, 39, 50, 50, 50}
	if logN == 16 {
		logQ[0], logQ[1] = 54, 38
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            logQ,
		LogP:            []int{50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	return params
}

func fastDFTBenchmarkMatrix(t testing.TB, params ckks.Parameters) Matrix {
	t.Helper()
	matrix, err := NewMatrixFromLiteral(params, MatrixLiteral{
		Type:         HomomorphicEncode,
		LogSlots:     4,
		LevelQ:       3,
		LevelP:       0,
		Levels:       []int{1, 1},
		Format:       Standard,
		LogBSGSRatio: 1,
	}, ckks.NewEncoder(params))
	require.NoError(t, err)
	return matrix
}

func BenchmarkFastDFTRepresentative(b *testing.B) {
	for _, logN := range []int{13, 16} {
		params := fastDFTBenchmarkParameters(b, logN)
		matrices := fastDFTBenchmarkMatrix(b, params)
		standardEval := standardDFTEvaluator(b, params, matrices)
		fastEval := NewFastEvaluator(params)
		standardInput := fastDFTInput(params, matrices.LevelQ)
		standardInput.LogDimensions = ring.Dimensions{Cols: matrices.LogSlots}
		fastInput := fastDFTMontgomeryCopy(params, standardInput)
		fastOut := ckks.NewCiphertext(params, 1, matrices.LevelQ)
		standardOut := ckks.NewCiphertext(params, 1, matrices.LevelQ)

		b.Run(fmt.Sprintf("LogN%d/Fast", logN), func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				fastOut.Resize(1, matrices.LevelQ)
				setFastOutputDomain(fastOut, fastInput)
				if err := fastEval.CoeffsToSlots(fastInput, matrices, fastOut, nil); err != nil {
					b.Fatal(err)
				}
			}
		})
		b.Run(fmt.Sprintf("LogN%d/Standard", logN), func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				standardOut.Resize(1, matrices.LevelQ)
				if err := standardEval.CoeffsToSlots(standardInput, matrices, standardOut, nil); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
