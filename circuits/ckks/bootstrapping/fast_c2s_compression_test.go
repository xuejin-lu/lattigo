package bootstrapping

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	ltcommon "github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"

	ckkslintrans "github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
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

	assertCheckpoint := func(name string, level int, scale rlwe.Scale, rowHash string) {
		t.Helper()
		require.Equal(t, level, state.Level(), "%s Level", name)
		require.Equal(t, 1, state.Degree(), "%s Degree", name)
		require.True(t, scale.Equal(state.Scale), "%s Scale: want %s, got %s", name, scale.Value.Text('e', 39), state.Scale.Value.Text('e', 39))
		require.Equal(t, rowHash, acceptedC2SRowHash(state, 3), "%s q0/q1/q2 rows", name)
	}

	wantScale := rlwe.NewScale(1 << 50)
	assertCheckpoint("group_0_input", 16, wantScale, "e4f6962ee08e26b1dc03a651c8a4975d4d68f6b364fc148b2a6f59d5129b46e4")
	checkpointHashes := [4][3]string{
		{"c7f05f9273347df4df086c850f619818d7402e6bbf7d8069069e20ad138a90e1", "c57cf3e340e158e701e6e56539e360d6eb860b4ca9d1450b78b286b3be54f22a", "2dc4d4242018adb06afd3cf779065b5db5597a721d2aa5ad3300b992efb46c4d"},
		{"2111e0e6bf829272ebab368bbf4f3de65a72fae62aebe74cd15eede4b57f8422", "5d613d674f2be122f6ae0b7227bf743382a44a53d6e3e11a27c7e75abfb1adb8", "f1a463babb9ef8a5c08e9df77d786f122d66c456153d8f6516fa8af290d48e43"},
		{"7f7e0ac133652266ca788467417f5e9d46b43b09e697e02b9e95620a7284f5ea", "8f49746a616a7b9fdbfb05fcdd0406b9f074123813992e4e43295e26356a3171", ""},
		{"faa566f804664ed30cd748c8998bc3bde2409c8d646bc21d2b6c7bc5e9adc004", "694627c930f6fd8aca8c1f80eefebf421b96203d124a3114924f8769a1548e48", ""},
	}
	matrixIndex := 0
	for group, factors := range eval.C2SDFTMatrix.Levels {
		for range factors {
			matrix := eval.C2SDFTMatrix.Matrices[matrixIndex]
			wantScale = wantScale.Mul(matrix.Scale)
			require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransform(state, ltcommon.LinearTransformation(matrix), state))
			matrixIndex++
		}
		levelBeforeRescale := state.Level()
		assertCheckpoint("group_raw", levelBeforeRescale, wantScale, checkpointHashes[group][0])

		wantScale = wantScale.Div(rlwe.NewScale(params.BootstrappingParameters.Q()[levelBeforeRescale]))
		require.NoError(t, eval.FastCKKS.Rescale(state, state))
		assertCheckpoint("group_post_rescale", levelBeforeRescale-1, wantScale, checkpointHashes[group][1])

		if exponent := eval.C2SRestorePlan[group]; exponent > 0 {
			factor := new(big.Int).Lsh(big.NewInt(1), uint(exponent))
			require.NoError(t, eval.FastCKKS.MulIntegerMaintained(state, factor, state))
			state.Scale = state.Scale.Mul(rlwe.NewScale(factor))
			wantScale = wantScale.Mul(rlwe.NewScale(factor))
			assertCheckpoint("group_post_restore", levelBeforeRescale-1, wantScale, checkpointHashes[group][2])
		}
	}
	require.Equal(t, len(eval.C2SDFTMatrix.Matrices), matrixIndex)
}

func TestFastModUpToLegacyC2SQ3AuthorityTransition(t *testing.T) {
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
	q012Before := make([][][]uint64, len(state.Value))
	for component := range state.Value {
		q3Before[component] = append([]uint64(nil), state.Value[component].Coeffs[3]...)
		q012Before[component] = make([][]uint64, 3)
		for row := 0; row < 3; row++ {
			q012Before[component][row] = append([]uint64(nil), state.Value[component].Coeffs[row]...)
		}
	}

	// QPREFIX-IMPL-006 has not migrated LinearTransform: this first legacy
	// C2S producer consumes q012 only. q3 remains physically populated but its
	// ModUp authority must no longer be assumed after this operation.
	firstMatrix := eval.C2SDFTMatrix.Matrices[0]
	require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransform(state, ltcommon.LinearTransformation(firstMatrix), state))
	q012Changed := false
	for component := range state.Value {
		require.Equal(t, q3Before[component], state.Value[component].Coeffs[3], "legacy C2S must not claim to preserve q3 authority")
		for row := 0; row < 3; row++ {
			if !equalUint64Slices(q012Before[component][row], state.Value[component].Coeffs[row]) {
				q012Changed = true
			}
		}
	}
	require.True(t, q012Changed, "the C2S producer must exercise its legacy q012 path")
}

func equalUint64Slices(a, b []uint64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestFastLogN13C2SProductionRescaleIgnoresPoisonedQ3(t *testing.T) {
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

	matrixIndex := 0
	for range eval.C2SDFTMatrix.Levels[0] {
		matrix := eval.C2SDFTMatrix.Matrices[matrixIndex]
		require.NoError(t, eval.DFTEvaluator.FastEvaluator().LinearTransform(state, ltcommon.LinearTransformation(matrix), state))
		matrixIndex++
	}
	poisoned := state.CopyNew()
	q3 := params.BootstrappingParameters.Q()[3]
	for component := range poisoned.Value {
		require.Len(t, poisoned.Value[component].Coeffs, state.Level()+1)
		require.Len(t, poisoned.Value[component].Coeffs[3], params.BootstrappingParameters.N())
		for i, residue := range poisoned.Value[component].Coeffs[3] {
			poisoned.Value[component].Coeffs[3][i] = (residue + uint64(101+i+component)) % q3
		}
	}

	cleanOut := ckks.NewCiphertext(params.BootstrappingParameters, 1, state.Level()-1)
	poisonedOut := ckks.NewCiphertext(params.BootstrappingParameters, 1, state.Level()-1)
	require.NoError(t, eval.FastCKKS.Rescale(state, cleanOut))
	require.NoError(t, eval.FastCKKS.Rescale(poisoned, poisonedOut))
	require.Equal(t, state.Level()-1, cleanOut.Level())
	require.True(t, cleanOut.Scale.Equal(poisonedOut.Scale))
	for component := range cleanOut.Value {
		for row := 0; row < 3; row++ {
			require.Equal(t, cleanOut.Value[component].Coeffs[row], poisonedOut.Value[component].Coeffs[row], "production C2S Rescale must ignore q3 component=%d q%d", component, row)
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
