package bootstrapping

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	ltcommon "github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"

	ckkslintrans "github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
	"github.com/tuneinsight/lattigo/v6/utils"
)

func fastLogN13CompressionParameters(t testing.TB) (Parameters, ckks.Parameters) {
	t.Helper()
	residual, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{56, 39},
		LogDefaultScale: 45,
		Xs:              ring.Ternary{H: 192},
	})
	require.NoError(t, err)
	zero := 0
	logN, logSlots := 13, 12
	params, err := NewParametersFromLiteral(residual, ParametersLiteral{
		LogN:     &logN,
		LogSlots: &logSlots,
		LogP:     []int{61, 61, 61, 61, 61},
		SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{39}, {39}, {39}},
		CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{56}, {56}, {56}, {56}},
		EvalModLogScale:       pointyInt(60),
		Mod1Degree:            pointyInt(30),
		DoubleAngle:           pointyInt(3),
		K:                     pointyInt(16),
		LogMessageRatio:       pointyInt(10),
		Mod1InvDegree:         &zero,
		EphemeralSecretWeight: &zero,
	})
	require.NoError(t, err)
	params.CircuitOrder = ModUpThenEncode
	params.ResidualParameters = residual
	return params, residual
}

func pointyInt(value int) *int { return &value }

func TestFastLogN13C2SCompressionPreparation(t *testing.T) {
	params, _ := fastLogN13CompressionParameters(t)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.NoError(t, eval.ensureFastBootstrapCircuit())
	require.True(t, eval.C2SCompressionActive)
	require.Equal(t, []int{4, 2, 0, 0}, eval.C2SRestorePlan)

	originalParams := params
	_, _, original, _, err := buildBootstrapCircuitData(originalParams)
	require.NoError(t, err)
	require.Equal(t, original.MatrixLiteral, eval.C2SDFTMatrix.MatrixLiteral)
	require.Equal(t, original.Levels, eval.C2SDFTMatrix.Levels)
	require.Len(t, eval.C2SDFTMatrix.Matrices, 4)
	for i := range original.Matrices {
		got := eval.C2SDFTMatrix.Matrices[i]
		want := original.Matrices[i]
		require.Equal(t, want.LevelQ, got.LevelQ)
		require.Equal(t, want.LevelP, got.LevelP)
		require.Equal(t, want.LogDimensions, got.LogDimensions)
		require.Equal(t, want.LogBabyStepGiantStepRatio, got.LogBabyStepGiantStepRatio)
		require.Equal(t, want.N1, got.N1)
		require.ElementsMatch(t, utils.GetKeys(want.Vec), utils.GetKeys(got.Vec))
		factor := 1.0
		if i == 0 {
			factor = 16
		} else if i == 1 {
			factor = 4
		}
		require.True(t, got.Scale.Mul(rlwe.NewScale(factor)).Equal(want.Scale), "group %d scale", i)
	}
}

func TestFastLogN13C2SCompressionGuardRejectsChangedProfile(t *testing.T) {
	params, _ := fastLogN13CompressionParameters(t)
	params.CoeffsToSlotsParameters.LogSlots = 11
	_, _, matrix, _, err := buildBootstrapCircuitData(params)
	require.NoError(t, err)
	require.False(t, matchesFastLogN13C2SProfile(params, matrix))
	prepared, plan, active, err := prepareFastLogN13C2S(params, matrix)
	require.NoError(t, err)
	require.False(t, active)
	require.Nil(t, plan)
	require.Equal(t, matrix, prepared)
}

func TestFastLogN13C2SCompressionUsesSameDiagonalStructure(t *testing.T) {
	params, _ := fastLogN13CompressionParameters(t)
	_, _, original, _, err := buildBootstrapCircuitData(params)
	require.NoError(t, err)
	prepared, plan, active, err := prepareFastLogN13C2S(params, original)
	require.NoError(t, err)
	require.True(t, active)
	require.Equal(t, []int{4, 2, 0, 0}, plan)
	mathematical := original.MatrixLiteral.GenMatrices(params.BootstrappingParameters.LogN(), params.BootstrappingParameters.EncodingPrecision())
	for i := range prepared.Matrices {
		require.ElementsMatch(t, mathematical[i].DiagonalsIndexList(), utils.GetKeys(original.Matrices[i].Vec), "group %d diagonal set", i)
		require.IsType(t, ckkslintrans.LinearTransformation{}, prepared.Matrices[i])
	}
}

func TestFastLogN13C2SRescaleMatchesAcceptedQPrefixAudit002Checkpoints(t *testing.T) {
	params, residual := fastLogN13CompressionParameters(t)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.NoError(t, eval.ensureFastBootstrapCircuit())
	require.Equal(t, []int{4, 2, 0, 0}, eval.C2SRestorePlan)

	values := make([]complex128, residual.MaxSlots())
	for i := range values {
		realValue := float64((i%7)-3) / 16
		imagValue := float64((i%5)-2) / 32
		values[i] = complex(realValue, imagValue)
	}
	state := fastBootstrapEncodedCiphertextAtLevel(t, residual, 0, params.CoeffsToSlotsParameters.LogSlots, values)
	state, _, err = eval.ScaleDown(state)
	require.NoError(t, err)
	state, err = eval.ModUp(state)
	require.NoError(t, err)
	require.Equal(t, []uint64{72057594037616641, 549755731969, 549756026881, 549755486209}, params.BootstrappingParameters.Q()[:4])

	assertCheckpoint := func(name string, level int, scale rlwe.Scale, rowHash string) []*big.Int {
		t.Helper()
		require.Equal(t, level, state.Level(), "%s Level", name)
		require.Equal(t, 1, state.Degree(), "%s Degree", name)
		require.True(t, scale.Equal(state.Scale), "%s Scale: want %s, got %s", name, scale.Value.Text('e', 39), state.Scale.Value.Text('e', 39))
		rows, err := fastckks.QPrefixWidth(state.Level())
		require.NoError(t, err)
		require.Equal(t, 4, rows, "%s authoritative width", name)
		gotHash := acceptedC2SRowHash(state, rows)
		if rowHash != "" {
			require.Equal(t, rowHash, gotHash, "%s q0/q1/q2/q3 rows", name)
		}
		bounds := exactCenteredQPrefixBounds(t, params.BootstrappingParameters, state, rows)
		product, err := fastckks.QPrefixProduct(params.BootstrappingParameters.Q(), state.Level())
		require.NoError(t, err)
		require.NoError(t, fastckks.CheckQPrefixCapacity(state.Level(), params.BootstrappingParameters.Q(), bounds), "%s strict capacity", name)
		rhos := make([]string, len(bounds))
		for component, bound := range bounds {
			numerator := new(big.Int).Lsh(new(big.Int).Set(bound), 1)
			rhos[component] = new(big.Rat).SetFrac(numerator, product).RatString()
		}
		t.Logf("%s Level=%d rows=%d Scale=%s Degree=%d q0123_product=%s max_abs=%v rho_2B_over_SQ=%v q0123_hash=%s", name, state.Level(), rows, state.Scale.Value.Text('e', 8), state.Degree(), product, bigIntStrings(bounds), rhos, gotHash)
		return bounds
	}

	wantScale := rlwe.NewScale(1 << 50)
	assertCheckpoint("group_0_input", 16, wantScale, "0919ab7a416e89ff4692192e1d23da5c4f05654d2b88d3ef8c5e42d6415be78e")
	checkpointHashes := [4][3]string{
		{"6de7db234649194726f6e4768eac94c4ba8ee575b2e3fda82d7cce19a6360b43", "9f9eed406067252265eee1859669276ed0b6fa5fc8d492ecc88692d06a6c34a7", "c2cb9363f5b1c3e8bcc19b1ff5bf3cbb98af8f6cf0e8d5313ad3f808f66563ae"},
		{"1bee43f083333d5a16b273f32e1a50b7d072a818cb93034088d01fe0d0481a30", "fa3cfebe842115b58c73e6cf59fffbae31585c2064a5b3e18681801dfe969668", "f2da9fe917c7528c19387283619d35f7c2c57f52d4f201636594c3fd20f012a5"},
		{"ada753851f80a75cf6f9a70512f2259e5ddd42ea7fb706127ab4a02a4582d89a", "192d1cc718a7c57cb66cf82033293a691797984bae1d410a3d28c8452c05f641", ""},
		{"80050cea635d89188101c61c33e64b1ac62d515efe557b31a330021723eef867", "03b608453b90b2a5b893d3200ffff4ba8995a37d34f1fb527e51460a5e9702f1", ""},
	}
	matrixIndex := 0
	for group, factors := range eval.C2SDFTMatrix.Levels {
		rows, err := fastckks.QPrefixWidth(state.Level())
		require.NoError(t, err)
		for range factors {
			matrix := eval.C2SDFTMatrix.Matrices[matrixIndex]
			wantScale = wantScale.Mul(matrix.Scale)
			require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransformQPrefixRows(state, ltcommon.LinearTransformation(matrix), rows, state))
			matrixIndex++
		}
		levelBeforeRescale := state.Level()
		rawBounds := assertCheckpoint("group_raw", levelBeforeRescale, wantScale, checkpointHashes[group][0])

		divisor := new(big.Int).SetUint64(params.BootstrappingParameters.Q()[levelBeforeRescale])
		wantScale = wantScale.Div(rlwe.NewScale(divisor))
		require.NoError(t, eval.FastCKKS.RescaleQPrefixRows(state, rows, state))
		targetRows, err := fastckks.QPrefixWidth(state.Level())
		require.NoError(t, err)
		rows = min(rows, targetRows)
		postBounds := assertCheckpoint("group_post_rescale", levelBeforeRescale-1, wantScale, checkpointHashes[group][1])
		rounding := new(big.Int).Rsh(new(big.Int).Sub(new(big.Int).Set(divisor), big.NewInt(1)), 1)
		for component, bound := range rawBounds {
			predicted := new(big.Int).Add(new(big.Int).Set(bound), rounding)
			predicted.Div(predicted, divisor)
			require.LessOrEqual(t, postBounds[component].Cmp(predicted), 0, "group %d component %d Rescale recurrence", group, component)
		}

		if exponent := eval.C2SRestorePlan[group]; exponent > 0 {
			factor := new(big.Int).Lsh(big.NewInt(1), uint(exponent))
			require.NoError(t, eval.FastCKKS.MulIntegerQPrefixRows(state, factor, rows, state))
			state.Scale = state.Scale.Mul(rlwe.NewScale(factor))
			wantScale = wantScale.Mul(rlwe.NewScale(factor))
			restoreBounds := assertCheckpoint("group_post_restore", levelBeforeRescale-1, wantScale, checkpointHashes[group][2])
			for component, bound := range postBounds {
				predicted := new(big.Int).Mul(new(big.Int).Set(bound), factor)
				require.LessOrEqual(t, restoreBounds[component].Cmp(predicted), 0, "group %d component %d restore recurrence", group, component)
			}
		}
	}
	require.Equal(t, len(eval.C2SDFTMatrix.Matrices), matrixIndex)
}

func TestFastModUpToProductionC2SPreservesQ3Authority(t *testing.T) {
	params, residual := fastLogN13CompressionParameters(t)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.NoError(t, eval.ensureFastBootstrapCircuit())
	values := make([]complex128, residual.MaxSlots())
	for i := range values {
		values[i] = complex(float64((i%7)-3)/16, float64((i%5)-2)/32)
	}
	state := fastBootstrapEncodedCiphertextAtLevel(t, residual, 0, params.CoeffsToSlotsParameters.LogSlots, values)
	state, _, err = eval.ScaleDown(state)
	require.NoError(t, err)
	state, err = eval.ModUp(state)
	require.NoError(t, err)
	require.Len(t, state.Value[0].Coeffs[3], params.BootstrappingParameters.N(), "ModUp establishes q3 authority")

	q3Before := make([][]uint64, len(state.Value))
	for component := range state.Value {
		q3Before[component] = append([]uint64(nil), state.Value[component].Coeffs[3]...)
	}

	// Production C2S must preserve ModUp's complete q0123 authority at the
	// first LinearTransform boundary.
	firstMatrix := eval.C2SDFTMatrix.Matrices[0]
	require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransformQPrefixRows(state, ltcommon.LinearTransformation(firstMatrix), 4, state))
	require.NotEqual(t, q3Before[0], state.Value[0].Coeffs[3], "C2S must transform c0 q3")
}

func TestFastLogN13C2SExplicitRowsConsumeQ3(t *testing.T) {
	params, residual := fastLogN13CompressionParameters(t)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.NoError(t, eval.ensureFastBootstrapCircuit())

	values := make([]complex128, residual.MaxSlots())
	for i := range values {
		values[i] = complex(float64((i%7)-3)/16, float64((i%5)-2)/32)
	}
	state := fastBootstrapEncodedCiphertextAtLevel(t, residual, 0, params.CoeffsToSlotsParameters.LogSlots, values)
	state, _, err = eval.ScaleDown(state)
	require.NoError(t, err)
	state, err = eval.ModUp(state)
	require.NoError(t, err)

	clean := state.CopyNew()
	poisoned := state.CopyNew()
	q3 := params.BootstrappingParameters.Q()[3]
	for component := range poisoned.Value {
		for i, residue := range poisoned.Value[component].Coeffs[3] {
			poisoned.Value[component].Coeffs[3][i] = (residue + uint64(101+i+component)) % q3
		}
	}
	firstMatrix := ltcommon.LinearTransformation(eval.C2SDFTMatrix.Matrices[0])
	require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransformQPrefixRows(clean, firstMatrix, 4, clean))
	require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransformQPrefixRows(poisoned, firstMatrix, 4, poisoned))
	for component := range clean.Value {
		for row := 0; row < 3; row++ {
			require.Equal(t, clean.Value[component].Coeffs[row], poisoned.Value[component].Coeffs[row], "q3 poison must not affect q%d component=%d", row, component)
		}
		require.NotEqual(t, clean.Value[component].Coeffs[3], poisoned.Value[component].Coeffs[3], "explicit production C2S must consume q3 component=%d", component)
	}
}

func TestFastLogN13S2CEvalModQPrefixBoundaryConsumesQ3(t *testing.T) {
	params, residual := fastLogN13CompressionParameters(t)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.NoError(t, eval.ensureFastBootstrapCircuit())
	values := make([]complex128, residual.MaxSlots())
	for i := range values {
		values[i] = complex(float64((i%7)-3)/16, float64((i%5)-2)/32)
	}
	input := fastBootstrapEncodedCiphertextAtLevel(t, residual, 0, params.CoeffsToSlotsParameters.LogSlots, values)
	input, _, err = eval.ScaleDown(input)
	require.NoError(t, err)
	input, err = eval.ModUp(input)
	require.NoError(t, err)
	ctReal, ctImag, err := eval.DFTEvaluator.CoeffsToSlotsNewWithRestorePlan(input, eval.C2SDFTMatrix, eval.C2SRestorePlan)
	require.NoError(t, err)
	ctReal, err = eval.EvalMod(ctReal)
	require.NoError(t, err)
	ctImag, err = eval.EvalMod(ctImag)
	require.NoError(t, err)
	require.GreaterOrEqual(t, ctReal.Level(), 3, "EvalMod fixture must have a physical q3 row")
	rows, err := fastckks.QPrefixWidth(ctReal.Level())
	require.NoError(t, err)
	require.Equal(t, 4, rows, "P93 EvalMod output must provide q0123 to production S2C")

	poison := func(ct *rlwe.Ciphertext) *rlwe.Ciphertext {
		out := ct.CopyNew()
		q3 := params.BootstrappingParameters.Q()[3]
		for component := range out.Value {
			for i, residue := range out.Value[component].Coeffs[3] {
				out.Value[component].Coeffs[3][i] = (residue + uint64(101+i+component)) % q3
			}
		}
		return out
	}
	cleanPrefix, poisonedPrefix := ctReal.CopyNew(), poison(ctReal)
	firstMatrix := ltcommon.LinearTransformation(eval.S2CDFTMatrix.Matrices[0])
	require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransformQPrefixRows(cleanPrefix, firstMatrix, rows, cleanPrefix))
	require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransformQPrefixRows(poisonedPrefix, firstMatrix, rows, poisonedPrefix))
	for component := range cleanPrefix.Value {
		for row := 0; row < 3; row++ {
			require.Equal(t, cleanPrefix.Value[component].Coeffs[row], poisonedPrefix.Value[component].Coeffs[row], "S2C q3 poison must not affect q%d component=%d", row, component)
		}
		require.NotEqual(t, cleanPrefix.Value[component].Coeffs[3], poisonedPrefix.Value[component].Coeffs[3], "the production S2C transform must consume q3 component=%d", component)
	}
	cleanOut, err := eval.SlotsToCoeffs(ctReal, ctImag)
	require.NoError(t, err)
	outputRows, err := fastckks.QPrefixWidth(cleanOut.Level())
	require.NoError(t, err)
	for component := range cleanOut.Value {
		require.GreaterOrEqual(t, len(cleanOut.Value[component].Coeffs), cleanOut.Level()+1)
		for row := 0; row < outputRows; row++ {
			require.Len(t, cleanOut.Value[component].Coeffs[row], params.BootstrappingParameters.N())
		}
		for row := outputRows; row <= cleanOut.Level(); row++ {
			require.Empty(t, cleanOut.Value[component].Coeffs[row], "S2C may contract q3 only when the logical Level contracts")
		}
	}
}

func acceptedC2SRowHash(ct *rlwe.Ciphertext, rows int) string {
	hash := sha256.New()
	var word [8]byte
	for component := range ct.Value {
		binary.LittleEndian.PutUint64(word[:], uint64(component))
		_, _ = hash.Write(word[:])
		for row := 0; row < rows; row++ {
			binary.LittleEndian.PutUint64(word[:], uint64(row))
			_, _ = hash.Write(word[:])
			for _, residue := range ct.Value[component].Coeffs[row] {
				binary.LittleEndian.PutUint64(word[:], residue)
				_, _ = hash.Write(word[:])
			}
		}
	}
	return hex.EncodeToString(hash.Sum(nil))
}

func exactCenteredQPrefixBounds(t testing.TB, params ckks.Parameters, ct *rlwe.Ciphertext, rows int) []*big.Int {
	t.Helper()
	require.GreaterOrEqual(t, rows, 1)
	require.LessOrEqual(t, rows, ct.Level()+1)
	require.LessOrEqual(t, rows, fastckks.MaxQPrefixWidth)
	q := params.Q()
	require.GreaterOrEqual(t, len(q), rows)
	products := make([]*big.Int, rows)
	inverses := make([]*big.Int, rows)
	prefix := big.NewInt(1)
	for row := 0; row < rows; row++ {
		products[row] = new(big.Int).Set(prefix)
		modulus := new(big.Int).SetUint64(q[row])
		prefixMod := new(big.Int).Mod(new(big.Int).Set(prefix), modulus)
		inverses[row] = new(big.Int).ModInverse(prefixMod, modulus)
		require.NotNil(t, inverses[row], "Q-prefix CRT inverse q%d", row)
		prefix.Mul(prefix, modulus)
	}
	half := new(big.Int).Rsh(new(big.Int).Set(prefix), 1)
	maxima := make([]*big.Int, len(ct.Value))
	ringQ := params.RingQ()
	for component := range ct.Value {
		coeffRows := make([][]uint64, rows)
		for row := 0; row < rows; row++ {
			coeffRows[row] = make([]uint64, params.N())
			copy(coeffRows[row], ct.Value[component].Coeffs[row])
			if ct.IsMontgomery {
				ringQ.SubRings[row].IMForm(coeffRows[row], coeffRows[row])
			}
			if ct.IsNTT {
				ringQ.SubRings[row].INTT(coeffRows[row], coeffRows[row])
			}
		}
		maxima[component] = new(big.Int)
		for coefficient := 0; coefficient < params.N(); coefficient++ {
			value := new(big.Int).SetUint64(coeffRows[0][coefficient])
			for row := 1; row < rows; row++ {
				modulus := new(big.Int).SetUint64(q[row])
				residue := new(big.Int).SetUint64(coeffRows[row][coefficient])
				current := new(big.Int).Mod(new(big.Int).Set(value), modulus)
				delta := new(big.Int).Sub(residue, current)
				delta.Mod(delta, modulus)
				factor := new(big.Int).Mul(delta, inverses[row])
				factor.Mod(factor, modulus)
				value.Add(value, new(big.Int).Mul(products[row], factor))
			}
			if value.Cmp(half) > 0 {
				value.Sub(value, prefix)
			}
			value.Abs(value)
			if value.Cmp(maxima[component]) > 0 {
				maxima[component].Set(value)
			}
		}
	}
	return maxima
}

func bigIntStrings(values []*big.Int) []string {
	strings := make([]string, len(values))
	for i, value := range values {
		strings[i] = value.String()
	}
	return strings
}

func BenchmarkFastLogN13P93C2SLinearTransformQPrefix(b *testing.B) {
	params, residual := fastLogN13CompressionParameters(b)
	eval, err := NewFastEvaluator(params)
	if err != nil {
		b.Fatal(err)
	}
	if err = eval.ensureFastBootstrapCircuit(); err != nil {
		b.Fatal(err)
	}
	values := make([]complex128, residual.MaxSlots())
	for i := range values {
		values[i] = complex(float64((i%7)-3)/16, float64((i%5)-2)/32)
	}
	input := fastBootstrapEncodedCiphertextAtLevel(b, residual, 0, params.CoeffsToSlotsParameters.LogSlots, values)
	input, _, err = eval.ScaleDown(input)
	if err != nil {
		b.Fatal(err)
	}
	input, err = eval.ModUp(input)
	if err != nil {
		b.Fatal(err)
	}

	seen := map[string]bool{}
	for matrixIndex, matrix := range eval.C2SDFTMatrix.Matrices {
		path := "direct"
		if matrix.N1 != 0 {
			path = "bsgs"
		}
		if seen[path] {
			continue
		}
		seen[path] = true
		lt := ltcommon.LinearTransformation(matrix)
		b.Run(fmt.Sprintf("group-matrix-%d/%s/legacy-q012", matrixIndex, path), func(b *testing.B) {
			out := input.CopyNew()
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				if err := eval.DFTEvaluator.FastEvaluator().LinearTransform(input, lt, out); err != nil {
					b.Fatal(err)
				}
			}
		})
		b.Run(fmt.Sprintf("group-matrix-%d/%s/explicit-q0123", matrixIndex, path), func(b *testing.B) {
			out := input.CopyNew()
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				if err := eval.DFTEvaluator.FastEvaluator().LinearTransformQPrefixRows(input, lt, 4, out); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
