package dft

import (
	"fmt"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	ckkslintrans "github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	ltcommon "github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

func fastDFTTestParameters(t testing.TB) ckks.Parameters {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		LogQ:            []int{56, 39, 39, 39, 50},
		LogP:            []int{50},
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	return params
}

func fastDFTMatrix(t testing.TB, params ckks.Parameters, typ Type, format Format) Matrix {
	t.Helper()
	return fastDFTMatrixAtLevel(t, params, typ, format, 4, []int{1, 1})
}

func fastDFTMatrixAtLevel(t testing.TB, params ckks.Parameters, typ Type, format Format, level int, groups []int) Matrix {
	t.Helper()
	literal := MatrixLiteral{
		Type:         typ,
		LogSlots:     2,
		LevelQ:       level,
		LevelP:       0,
		Levels:       groups,
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

func requireFastDFTMatchesStandard(t *testing.T, params ckks.Parameters, fast *rlwe.Ciphertext, standard *rlwe.Ciphertext, rows int) {
	t.Helper()
	got := fastDFTNonMontgomeryCopy(params, fast)
	require.Equal(t, standard.Level(), got.Level())
	require.Equal(t, standard.Scale, got.Scale)
	require.Equal(t, standard.LogDimensions, got.LogDimensions)
	for d := range standard.Value {
		require.Equal(t, standard.Value[d].Coeffs[:rows], got.Value[d].Coeffs[:rows], "component %d", d)
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
	requireFastDFTMatchesStandard(t, params, fastReal, standardReal, min(4, standardReal.Level()+1))
	require.Equal(t, cts.LevelQ-len(cts.Levels)*params.LevelsConsumedPerRescaling(), fastReal.Level())

	fastDecoded, err := fastEval.SlotsToCoeffsNew(fastInput, nil, stc)
	require.NoError(t, err)
	standardDecoded, err := standardEval.SlotsToCoeffsNew(standardInput, nil, stc)
	require.NoError(t, err)
	s2cRows := min(fastckks.MaintainedLimbCount(&params, fastInput.Level()), standardDecoded.Level()+1)
	requireFastDFTMatchesStandard(t, params, fastDecoded, standardDecoded, s2cRows)
}

func TestFastDFTCoeffsToSlotsRestorePlan(t *testing.T) {
	params := fastDFTTestParameters(t)
	original := fastDFTMatrix(t, params, HomomorphicEncode, Standard)
	compressed := original
	compressed.Matrices = append([]ckkslintrans.LinearTransformation(nil), original.Matrices...)
	mathematical := original.MatrixLiteral.GenMatrices(params.LogN(), params.EncodingPrecision())
	for i, k := range []int{2, 0} {
		if k == 0 {
			continue
		}
		matrix := original.Matrices[i]
		compressed.Matrices[i] = ckkslintrans.NewTransformation(params, ckkslintrans.Parameters{
			DiagonalsIndexList:        mathematical[i].DiagonalsIndexList(),
			LevelQ:                    matrix.LevelQ,
			LevelP:                    matrix.LevelP,
			Scale:                     matrix.Scale.Div(rlwe.NewScale(uint64(1) << uint(k))),
			LogDimensions:             matrix.LogDimensions,
			LogBabyStepGiantStepRatio: matrix.LogBabyStepGiantStepRatio,
		})
		require.NoError(t, ckkslintrans.Encode(ckks.NewEncoder(params), mathematical[i], compressed.Matrices[i]))
	}

	standardInput := fastDFTInput(params, original.LevelQ)
	fastInput := fastDFTMontgomeryCopy(params, standardInput)
	standardEval := standardDFTEvaluator(t, params, compressed)
	standard := standardInput.CopyNew()
	matrixIdx := 0
	for groupIdx, factors := range compressed.Levels {
		for range factors {
			out := ckks.NewCiphertext(params, 1, standard.Level())
			require.NoError(t, standardEval.LTEvaluator.Evaluate(standard, compressed.Matrices[matrixIdx], out))
			standard = out
			matrixIdx++
		}
		require.NoError(t, standardEval.Rescale(standard, standard))
		if k := []int{2, 0}[groupIdx]; k != 0 {
			require.NoError(t, standardEval.Mul(standard, uint64(1)<<uint(k), standard))
			standard.Scale = standard.Scale.Mul(rlwe.NewScale(uint64(1) << uint(k)))
		}
	}
	fast, _, err := NewFastEvaluator(params).CoeffsToSlotsNewWithRestorePlan(fastInput, compressed, []int{2, 0})
	require.NoError(t, err)
	requireFastDFTMatchesStandard(t, params, fast, standard, min(4, standard.Level()+1))
}

func TestFastDFTS2CQPrefixCapabilityMatchesStandard(t *testing.T) {
	params := fastDFTTestParameters(t)
	matrices := fastDFTMatrixAtLevel(t, params, HomomorphicDecode, Standard, 3, []int{1, 1})
	fastEval := NewFastEvaluator(params)
	standardEval := standardDFTEvaluator(t, params, matrices)
	standardInput := fastDFTInput(params, 3)
	fastInput := fastDFTMontgomeryCopy(params, standardInput)
	standardImag := fastDFTInput(params, 3)
	fastImag := fastDFTMontgomeryCopy(params, standardImag)

	// Compare the first factor at all four rows while q3 is authoritative.
	fastFirst := ckks.NewCiphertext(params, 1, 3)
	fastFirst.IsNTT, fastFirst.IsMontgomery = true, true
	require.NoError(t, fastEval.FastEvaluator().LinearTransformQPrefixRows(fastInput, ltcommon.LinearTransformation(matrices.Matrices[0]), 4, fastFirst))
	standardFirst := ckks.NewCiphertext(params, 1, 3)
	require.NoError(t, standardEval.LTEvaluator.Evaluate(standardInput, matrices.Matrices[0], standardFirst))
	gotFirst := fastDFTNonMontgomeryCopy(params, fastFirst)
	for d := range standardFirst.Value {
		require.Equal(t, standardFirst.Value[d].Coeffs[:4], gotFirst.Value[d].Coeffs[:4], "first factor component=%d", d)
	}

	fastOut, err := fastEval.SlotsToCoeffsNewQPrefixRows(fastInput, fastImag, matrices, 4)
	require.NoError(t, err)
	rows, err := fastckks.QPrefixWidth(fastOut.Level())
	require.NoError(t, err)
	standardOut, err := standardEval.SlotsToCoeffsNew(standardInput, standardImag, matrices)
	require.NoError(t, err)
	require.Equal(t, 2, rows, "3->2->1 natural contraction")
	requireFastDFTMatchesStandard(t, params, fastOut, standardOut, 2)
}

func TestFastDFTS2CProductionDoesNotPromoteEvalModQ3(t *testing.T) {
	params := fastDFTTestParameters(t)
	matrices := fastDFTMatrixAtLevel(t, params, HomomorphicDecode, Standard, 3, []int{1, 1})
	fastEval := NewFastEvaluator(params)
	cleanInput := fastDFTMontgomeryCopy(params, fastDFTInput(params, 3))
	poisonedInput := cleanInput.CopyNew()
	q3 := params.RingQ().SubRings[3].Modulus
	for d := range poisonedInput.Value {
		for i, value := range poisonedInput.Value[d].Coeffs[3] {
			poisonedInput.Value[d].Coeffs[3][i] = (value + uint64(101+i+d)) % q3
		}
	}

	clean, err := fastEval.SlotsToCoeffsNew(cleanInput, nil, matrices)
	require.NoError(t, err)
	poisoned, err := fastEval.SlotsToCoeffsNew(poisonedInput, nil, matrices)
	require.NoError(t, err)
	legacyRows := fastckks.MaintainedLimbCount(&params, cleanInput.Level())
	require.Less(t, legacyRows, 4, "production S2C must not promote allocated q3")
	require.Equal(t, clean.Level(), poisoned.Level())
	require.Equal(t, clean.Scale, poisoned.Scale)
	for d := range clean.Value {
		for row := 0; row < min(legacyRows, clean.Level()+1); row++ {
			require.Equal(t, clean.Value[d].Coeffs[row], poisoned.Value[d].Coeffs[row], "S2C consumed poisoned non-authoritative q3 component=%d q%d", d, row)
		}
	}
}

func TestFastDFTS2CQPrefixRescaleContractsThroughLevelZero(t *testing.T) {
	params := fastDFTTestParameters(t)
	levels := []int{1, 1, 1}
	matrices, err := NewMatrixFromLiteral(params, MatrixLiteral{
		Type:         HomomorphicDecode,
		LogSlots:     2,
		LevelQ:       3,
		LevelP:       0,
		Levels:       levels,
		Format:       Standard,
		LogBSGSRatio: 1,
	}, ckks.NewEncoder(params))
	require.NoError(t, err)
	fastEval := NewFastEvaluator(params)
	standardEval := standardDFTEvaluator(t, params, matrices)
	standardInput := fastDFTInput(params, 3)
	fastInput := fastDFTMontgomeryCopy(params, standardInput)
	fastOut, err := fastEval.SlotsToCoeffsNewQPrefixRows(fastInput, nil, matrices, 4)
	require.NoError(t, err)
	rows, err := fastckks.QPrefixWidth(fastOut.Level())
	require.NoError(t, err)
	standardOut, err := standardEval.SlotsToCoeffsNew(standardInput, nil, matrices)
	require.NoError(t, err)
	require.Equal(t, 0, fastOut.Level(), "three groups contract 3->2->1->0")
	require.Equal(t, 1, rows)
	requireFastDFTMatchesStandard(t, params, fastOut, standardOut, 1)
}

func TestFastDFTSplitRepackAndPoison(t *testing.T) {
	params := fastDFTTestParameters(t)
	matrix := fastDFTMatrix(t, params, HomomorphicEncode, RepackImagAsReal)
	fastEval := NewFastEvaluator(params)
	in := fastDFTMontgomeryCopy(params, fastDFTInput(params, matrix.LevelQ))
	poisoned := in.CopyNew()
	for d := range poisoned.Value {
		for limb := 4; limb <= poisoned.Level(); limb++ {
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
	require.Equal(t, fastReal.Value[0].Coeffs[:3], poisonedReal.Value[0].Coeffs[:3])
	require.Equal(t, fastReal.Value[1].Coeffs[:3], poisonedReal.Value[1].Coeffs[:3])
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
	b.Run("FastQPrefix", func(b *testing.B) {
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
