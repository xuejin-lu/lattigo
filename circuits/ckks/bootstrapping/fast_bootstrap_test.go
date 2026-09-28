package bootstrapping

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
	"github.com/tuneinsight/lattigo/v6/utils"
)

func fastBootstrapParameters(t testing.TB, logSlots int, factorTwo bool) (Parameters, ckks.Parameters) {
	t.Helper()
	logN := 4
	return fastBootstrapParametersAt(t, logN, logSlots, factorTwo)
}

func fastBootstrapParametersAt(t testing.TB, logN, logSlots int, factorTwo bool) (Parameters, ckks.Parameters) {
	return fastBootstrapParametersProfileAt(t, logN, logSlots, factorTwo, 4, 30, 2)
}

func fastBootstrapParametersProfileAt(t testing.TB, logN, logSlots int, factorTwo bool, k, mod1Degree, doubleAngle int) (Parameters, ckks.Parameters) {
	t.Helper()
	residualLogN := logN
	if factorTwo {
		residualLogN = logN - 1
	}
	residualLiteral := ckks.ParametersLiteral{LogN: residualLogN, LogQ: []int{55, 39}, LogDefaultScale: 30}
	if logN == 16 {
		residualLiteral.LogQ = []int{54, 38}
	}
	if factorTwo {
		residualLiteral.LogNthRoot = logN + 1
	}
	residual, err := ckks.NewParametersFromLiteral(residualLiteral)
	require.NoError(t, err)
	zero := 0
	factorization := [][]int{{39}}
	if logSlots > 1 {
		factorization = [][]int{{39}, {39}}
	}
	literal := ParametersLiteral{
		LogN:     utils.Pointy(logN),
		LogSlots: utils.Pointy(logSlots),
		SlotsToCoeffsFactorizationDepthAndLogScales: factorization,
		CoeffsToSlotsFactorizationDepthAndLogScales: factorization,
		EvalModLogScale:       utils.Pointy(45),
		Mod1Degree:            utils.Pointy(mod1Degree),
		DoubleAngle:           utils.Pointy(doubleAngle),
		K:                     utils.Pointy(k),
		LogMessageRatio:       utils.Pointy(2),
		Mod1InvDegree:         &zero,
		EphemeralSecretWeight: &zero,
	}
	params, err := NewParametersFromLiteral(residual, literal)
	require.NoError(t, err)
	return params, residual
}

func actualLogN16BootstrapParameters(t testing.TB) (Parameters, ckks.Parameters) {
	t.Helper()
	residual, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            16,
		LogQ:            []int{55, 39},
		LogDefaultScale: 45,
		Xs:              ring.Ternary{H: 192},
	})
	require.NoError(t, err)
	logN := 16
	logSlots := residual.LogMaxSlots()
	evalModScale := 60
	mod1Degree := 30
	doubleAngle := 3
	k := 16
	logMessageRatio := 10
	invDegree := 0
	zero := 0
	params, err := NewParametersFromLiteral(residual, ParametersLiteral{
		LogN:     &logN,
		LogP:     []int{61, 61, 61, 61, 61},
		Xs:       ring.Ternary{H: 192},
		LogSlots: &logSlots,
		CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{56}, {56}, {56}, {56}},
		SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{39}, {39}, {39}},
		EvalModLogScale:       &evalModScale,
		EphemeralSecretWeight: &zero,
		Mod1Type:              mod1.CosDiscrete,
		LogMessageRatio:       &logMessageRatio,
		K:                     &k,
		Mod1Degree:            &mod1Degree,
		DoubleAngle:           &doubleAngle,
		Mod1InvDegree:         &invDegree,
	})
	require.NoError(t, err)
	params.CircuitOrder = ModUpThenEncode
	return params, residual
}

func fastBootstrapPlainCiphertext(t testing.TB, params ckks.Parameters, level, logSlots int) *rlwe.Ciphertext {
	t.Helper()
	ct := ckks.NewCiphertext(params, 1, level)
	ct.IsNTT = true
	ct.IsMontgomery = false
	ct.IsBatched = true
	ct.Scale = params.DefaultScale()
	ct.LogDimensions = ring.Dimensions{Cols: logSlots}
	for d := range ct.Value {
		for limb := 0; limb <= level; limb++ {
			for i := range ct.Value[d].Coeffs[limb] {
				value := int64((i % 5) - 2)
				if d == 1 {
					value = 0
				}
				ct.Value[d].Coeffs[limb][i] = new(big.Int).Mod(big.NewInt(value), new(big.Int).SetUint64(params.RingQ().SubRings[limb].Modulus)).Uint64()
			}
		}
		params.RingQ().AtLevel(level).NTT(ct.Value[d], ct.Value[d])
	}
	return ct
}

func fastBootstrapEncodedCiphertext(t testing.TB, params ckks.Parameters, logSlots int, values []complex128) *rlwe.Ciphertext {
	return fastBootstrapEncodedCiphertextAtLevel(t, params, 0, logSlots, values)
}

func fastBootstrapEncodedCiphertextAtLevel(t testing.TB, params ckks.Parameters, level, logSlots int, values []complex128) *rlwe.Ciphertext {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	pt := ckks.NewPlaintext(params, level)
	pt.IsNTT = true
	pt.IsMontgomery = false
	pt.LogDimensions = ring.Dimensions{Cols: logSlots}
	require.NoError(t, encoder.Encode(values, pt))
	ct := ckks.NewCiphertext(params, 1, level)
	*ct.MetaData = *pt.MetaData
	ct.Value[0].Copy(pt.Value)
	ct.Value[1].Zero()
	ct.IsNTT = pt.IsNTT
	ct.IsMontgomery = pt.IsMontgomery
	return ct
}

func fastBootstrapDecode(t *testing.T, params ckks.Parameters, ct *rlwe.Ciphertext, values []complex128) {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	pt := ckks.NewPlaintext(params, ct.Level())
	*pt.MetaData = *ct.MetaData
	pt.Value.Copy(ct.Value[0])
	pt.IsNTT = ct.IsNTT
	pt.IsMontgomery = ct.IsMontgomery
	decoded := make([]complex128, len(values))
	require.NoError(t, encoder.Decode(pt, decoded))
	for i := range values {
		require.InDelta(t, real(values[i]), real(decoded[i]), 1e-2, "real slot %d", i)
		require.InDelta(t, imag(values[i]), imag(decoded[i]), 1e-2, "imag slot %d", i)
	}
}

func TestFastBootstrapCircuitInitializationAndInterface(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.Equal(t, 0, eval.MinimumInputLevel())
	require.Equal(t, residual.MaxLevel(), eval.OutputLevel())
	require.NotNil(t, eval)
	require.NoError(t, eval.ensureFastBootstrapCircuit())
	require.NotNil(t, eval.DFTEvaluator)
	require.NotEmpty(t, eval.C2SDFTMatrix.Matrices)
	require.NotEmpty(t, eval.S2CDFTMatrix.Matrices)
	require.Equal(t, mod1.CosDiscrete, eval.Mod1Parameters.Mod1Type)
	require.Nil(t, eval.Mod1Parameters.Mod1InvPoly)
}

func TestFastBootstrapPublicValidation(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	_, err = eval.Bootstrap(nil)
	require.Error(t, err)
	_, err = eval.BootstrapMany(nil)
	require.Error(t, err)

	ct := fastBootstrapPlainCiphertext(t, residual, 0, 2)
	ct.IsMontgomery = true
	_, err = eval.Bootstrap(ct)
	require.Error(t, err)
}

func TestFastBootstrapStageAEndToEndSmoke(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	ct := fastBootstrapPlainCiphertext(t, residual, 0, 2)
	source := ct.CopyNew()
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	require.Equal(t, residual.N(), out.N())
	require.Equal(t, residual.MaxLevel(), out.Level())
	require.True(t, out.IsNTT)
	require.False(t, out.IsMontgomery)
	require.True(t, out.Scale.Equal(residual.DefaultScale()))
	require.Equal(t, source.Value[0].Coeffs[0], ct.Value[0].Coeffs[0])
}

func TestFastBootstrapPreservesNonZeroMessage(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i}
	ct := fastBootstrapEncodedCiphertext(t, residual, 2, values)
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	fastBootstrapDecode(t, residual, out, values)
}

func TestFastBootstrapLevelOneInput(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i}
	ct := fastBootstrapEncodedCiphertextAtLevel(t, residual, 1, 2, values)
	source := ct.CopyNew()
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	require.Equal(t, residual.MaxLevel(), out.Level())
	require.Equal(t, source.Value[0].Coeffs[:2], ct.Value[0].Coeffs[:2])
	fastBootstrapDecode(t, residual, out, values)
}

func TestFastBootstrapFullSlotRealImagEndToEnd(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 3, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{
		0.03125 + 0.0625i,
		-0.0625 + 0.09375i,
		0.125 - 0.03125i,
		-0.15625 - 0.125i,
		0.1875 + 0.15625i,
		-0.21875 + 0.1875i,
		0.25 - 0.21875i,
		-0.28125 - 0.25i,
	}
	ct := fastBootstrapEncodedCiphertext(t, residual, 3, values)
	require.NoError(t, eval.ensureFastBootstrapCircuit())

	// Exercise the same public-input boundary as Bootstrap before the public
	// call and verify that the full-slot C2S branch produces both outputs.
	coreInput := ct.CopyNew()
	coreInput, _, err = eval.ScaleDown(coreInput)
	require.NoError(t, err)
	coreInput, err = eval.ModUp(coreInput)
	require.NoError(t, err)
	ctReal, ctImag, err := eval.DFTEvaluator.CoeffsToSlotsNew(coreInput, eval.C2SDFTMatrix)
	require.NoError(t, err)
	require.NotNil(t, ctReal)
	require.NotNil(t, ctImag)

	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	require.True(t, out.IsNTT)
	require.False(t, out.IsMontgomery)
	require.Equal(t, residual.MaxLevel(), out.Level())
	require.True(t, out.Scale.Equal(residual.DefaultScale()))
	fastBootstrapDecode(t, residual, out, values)
}

func TestFastBootstrapActualLogN16FullSlotProfile(t *testing.T) {
	params, residual := actualLogN16BootstrapParameters(t)
	params.ResidualParameters = residual
	require.Equal(t, 15, params.LogMaxSlots())
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := make([]complex128, residual.MaxSlots())
	for i := range values {
		values[i] = complex(float64((i%7)-3)/16, float64((i%5)-2)/32)
	}
	ct := fastBootstrapEncodedCiphertext(t, residual, params.LogMaxSlots(), values)
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	require.Equal(t, residual.MaxLevel(), out.Level())
	require.True(t, out.IsNTT)
	require.False(t, out.IsMontgomery)
}

func TestFastBootstrapSparseRepackedPath(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 1, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i}
	ct := fastBootstrapEncodedCiphertext(t, residual, 1, values)
	out, err := eval.Bootstrap(ct)
	require.NoError(t, err)
	fastBootstrapDecode(t, residual, out, values)
}

func TestFastBootstrapManyOddCount(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := make([][]complex128, 3)
	cts := make([]rlwe.Ciphertext, len(values))
	for i := range cts {
		values[i] = []complex128{complex(0.05*float64(i+1), 0.01), complex(-0.04*float64(i+1), -0.02), complex(0.03+0.01*float64(i), 0.05), complex(-0.02, -0.04-0.01*float64(i))}
		cts[i] = *fastBootstrapEncodedCiphertext(t, residual, 2, values[i])
	}
	outputs, err := eval.BootstrapMany(cts)
	require.NoError(t, err)
	require.Len(t, outputs, 3)
	for i, out := range outputs {
		require.Equal(t, 2, out.LogSlots())
		require.False(t, out.IsMontgomery)
		require.Equal(t, residual.MaxLevel(), out.Level())
		fastBootstrapDecode(t, residual, &out, values[i])
	}
	other := outputs[1].Value[0].Coeffs[0][0]
	outputs[0].Value[0].Coeffs[0][0]++
	require.Equal(t, other, outputs[1].Value[0].Coeffs[0][0], "BootstrapMany outputs must not alias")
}

func TestFastBootstrapManyStructuralMatrix(t *testing.T) {
	cases := []struct {
		name      string
		count     int
		factorTwo bool
		level     int
		logSlots  int
	}{
		{name: "count1-level0-n1eqn2", count: 1, level: 0, logSlots: 2},
		{name: "count2-level1-n1eqn2", count: 2, level: 1, logSlots: 2},
		{name: "count3-sparse-n1eqn2", count: 3, level: 0, logSlots: 1},
		{name: "count5-level1-n1eqn2", count: 5, level: 1, logSlots: 2},
		{name: "count3-level0-factor-two", count: 3, factorTwo: true, level: 0, logSlots: 2},
		{name: "count5-level1-factor-two", count: 5, factorTwo: true, level: 1, logSlots: 2},
		{name: "count3-full-slots", count: 3, level: 0, logSlots: 3},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var params Parameters
			var residual ckks.Parameters
			if tc.factorTwo {
				params, residual = fastBootstrapParametersAt(t, 5, tc.logSlots, true)
			} else {
				params, residual = fastBootstrapParameters(t, tc.logSlots, false)
			}
			params.ResidualParameters = residual
			eval, err := NewFastEvaluator(params)
			require.NoError(t, err)
			values := make([][]complex128, tc.count)
			inputs := make([]rlwe.Ciphertext, tc.count)
			original := make([]rlwe.Ciphertext, tc.count)
			for i := range inputs {
				values[i] = make([]complex128, 1<<tc.logSlots)
				for j := range values[i] {
					values[i][j] = complex(float64((i+j)%5-2)/32, float64((2*i+j)%7-3)/64)
				}
				input := fastBootstrapEncodedCiphertextAtLevel(t, residual, tc.level, tc.logSlots, values[i])
				inputs[i] = *input
				original[i] = *input.CopyNew()
			}

			outputs, err := eval.BootstrapMany(inputs)
			require.NoError(t, err)
			require.Len(t, outputs, tc.count)
			rows, err := fastckks.QPrefixWidth(residual.MaxLevel())
			require.NoError(t, err)
			for i := range outputs {
				out := &outputs[i]
				require.Equal(t, residual.N(), out.N())
				require.Equal(t, residual.MaxLevel(), out.Level())
				require.Equal(t, 1, out.Degree())
				require.True(t, out.IsNTT)
				require.False(t, out.IsMontgomery)
				require.True(t, out.Scale.Equal(residual.DefaultScale()))
				require.Equal(t, tc.logSlots, out.LogSlots())
				require.Len(t, out.Value, 2)
				for d := 0; d <= 1; d++ {
					require.Len(t, out.Value[d].Coeffs, rows)
					for row := 0; row < rows; row++ {
						require.Len(t, out.Value[d].Coeffs[row], residual.N())
					}
				}
				fastBootstrapDecode(t, residual, out, values[i])
				require.Equal(t, *original[i].MetaData, *inputs[i].MetaData, "input %d metadata mutated", i)
				for d := 0; d <= 1; d++ {
					require.Equal(t, original[i].Value[d].Coeffs, inputs[i].Value[d].Coeffs, "input %d component %d mutated", i, d)
				}
			}
			if len(outputs) > 1 {
				other := outputs[1].Value[0].Coeffs[0][0]
				outputs[0].Value[0].Coeffs[0][0]++
				require.Equal(t, other, outputs[1].Value[0].Coeffs[0][0], "BootstrapMany outputs must not alias")
			}
		})
	}
}

func TestFastBootstrapPublicBoundaryRejectsRowsBeyondContract(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i}
	base := fastBootstrapEncodedCiphertextAtLevel(t, residual, 1, 2, values)
	poisoned := base.CopyNew()
	for d := 0; d <= 1; d++ {
		for row := 2; row <= 3; row++ {
			poisoned.Value[d].Coeffs = append(poisoned.Value[d].Coeffs, make([]uint64, residual.N()))
			for i := range poisoned.Value[d].Coeffs[row] {
				poisoned.Value[d].Coeffs[row][i] = uint64(0x9e3779b9 + 31*d + 17*row + i)
			}
		}
	}
	baseBefore, poisonedBefore := base.CopyNew(), poisoned.CopyNew()
	_, err = eval.Bootstrap(base)
	require.NoError(t, err)
	_, err = eval.Bootstrap(poisoned)
	require.Error(t, err, "physical extra rows increase the RLWE logical Level and violate the public residual contract")
	for d := 0; d <= 1; d++ {
		require.Equal(t, baseBefore.Value[d].Coeffs, base.Value[d].Coeffs, "base input component %d mutated", d)
		require.Equal(t, poisonedBefore.Value[d].Coeffs, poisoned.Value[d].Coeffs, "poisoned input component %d mutated", d)
	}
}

func TestFastBootstrapFinalizationRejectsMalformedRowsTransactionally(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	ct := ckks.NewCiphertext(residual, 1, residual.MaxLevel())
	ct.IsNTT, ct.IsMontgomery = true, true
	ct.Value[0].Coeffs[1] = ct.Value[0].Coeffs[1][:residual.N()-1]
	before := ct.CopyNew()
	require.Error(t, eval.finalizeFastPublicCiphertext(ct))
	require.Equal(t, before.Value[0].Coeffs, ct.Value[0].Coeffs)
	require.Equal(t, before.Value[1].Coeffs, ct.Value[1].Coeffs)
	require.Equal(t, before.IsNTT, ct.IsNTT)
	require.Equal(t, before.IsMontgomery, ct.IsMontgomery)
}

func TestFastBootstrapManyMatchesExplicitQPrefixStageSequence(t *testing.T) {
	params, residual := fastBootstrapParametersAt(t, 5, 2, true)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i}
	input := fastBootstrapEncodedCiphertextAtLevel(t, residual, 1, 2, values)
	public, err := eval.Bootstrap(input.CopyNew())
	require.NoError(t, err)

	inputRows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)
	packed, ctxtN1, ctxtN2, err := eval.PackAndSwitchN1ToN2QPrefixRows([]rlwe.Ciphertext{*input.CopyNew()}, inputRows)
	require.NoError(t, err)
	staged := &packed[0]
	staged, _, err = eval.ScaleDown(staged)
	require.NoError(t, err)
	staged, err = eval.ModUp(staged)
	require.NoError(t, err)
	realPart, imagPart, err := eval.DFTEvaluator.CoeffsToSlotsNewWithRestorePlan(staged, eval.C2SDFTMatrix, eval.C2SRestorePlan)
	require.NoError(t, err)
	realPart, err = eval.EvalMod(realPart)
	require.NoError(t, err)
	if imagPart != nil {
		imagPart, err = eval.EvalMod(imagPart)
		require.NoError(t, err)
	}
	staged, err = eval.SlotsToCoeffs(realPart, imagPart)
	require.NoError(t, err)
	outputRows, err := fastckks.QPrefixWidth(staged.Level())
	require.NoError(t, err)
	unpacked, err := eval.UnpackAndSwitchN2ToN1QPrefixRows([]rlwe.Ciphertext{*staged}, ctxtN1, ctxtN2, outputRows)
	require.NoError(t, err)
	require.Len(t, unpacked, 1)
	require.NoError(t, eval.finalizeFastPublicCiphertext(&unpacked[0]))
	requireFastBootstrapPublicEqual(t, public, &unpacked[0], residual)
	fastBootstrapDecode(t, residual, &unpacked[0], values)
}

func requireFastBootstrapPublicEqual(t *testing.T, want, got *rlwe.Ciphertext, params ckks.Parameters) {
	t.Helper()
	require.Equal(t, params.N(), got.N())
	require.Equal(t, params.MaxLevel(), got.Level())
	require.Equal(t, want.Degree(), got.Degree())
	require.Equal(t, *want.MetaData, *got.MetaData)
	require.Equal(t, want.IsNTT, got.IsNTT)
	require.Equal(t, want.IsMontgomery, got.IsMontgomery)
	require.Equal(t, want.Value[0].Coeffs, got.Value[0].Coeffs)
	require.Equal(t, want.Value[1].Coeffs, got.Value[1].Coeffs)
}

func TestFastBootstrapOrdinaryDecryptorDecode(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i}
	input := fastBootstrapEncodedCiphertext(t, residual, 2, values)
	out, err := eval.Bootstrap(input)
	require.NoError(t, err)

	zeroSK := rlwe.NewSecretKey(residual)
	zeroSK.Value.Q.Zero()
	decryptor := rlwe.NewDecryptor(residual, zeroSK)
	plaintext := decryptor.DecryptNew(out)
	decoded := make([]complex128, len(values))
	require.NoError(t, ckks.NewEncoder(residual).Decode(plaintext, decoded))
	for i := range values {
		require.InDelta(t, real(values[i]), real(decoded[i]), 1e-2, "real slot %d", i)
		require.InDelta(t, imag(values[i]), imag(decoded[i]), 1e-2, "imag slot %d", i)
	}
}

func TestFastBootstrapPlanningMatchesStandard(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.NoError(t, fastEval.ensureFastBootstrapCircuit())

	sk := rlwe.NewKeyGenerator(params.BootstrappingParameters).GenSecretKeyNew()
	keys, _, err := params.GenEvaluationKeys(sk)
	require.NoError(t, err)
	standardEval, err := NewEvaluator(params, keys)
	require.NoError(t, err)

	require.Equal(t, standardEval.Mod1Parameters.LevelQ, fastEval.Mod1Parameters.LevelQ)
	require.Equal(t, standardEval.Mod1Parameters.QDiff, fastEval.Mod1Parameters.QDiff)
	for _, matrices := range [][2]dft.Matrix{{standardEval.C2SDFTMatrix, fastEval.C2SDFTMatrix}, {standardEval.S2CDFTMatrix, fastEval.S2CDFTMatrix}} {
		require.Equal(t, matrices[0].LevelQ, matrices[1].LevelQ)
		require.Equal(t, matrices[0].Levels, matrices[1].Levels)
		require.Equal(t, matrices[0].LogSlots, matrices[1].LogSlots)
		require.NotNil(t, matrices[0].Scaling)
		require.NotNil(t, matrices[1].Scaling)
		require.Zero(t, matrices[0].Scaling.Cmp(matrices[1].Scaling))
	}
	require.NotNil(t, standardEval.Parameters.CoeffsToSlotsParameters.Scaling)
	require.NotNil(t, fastEval.Parameters.CoeffsToSlotsParameters.Scaling)
	require.Zero(t, standardEval.Parameters.CoeffsToSlotsParameters.Scaling.Cmp(fastEval.Parameters.CoeffsToSlotsParameters.Scaling))
	require.NotNil(t, standardEval.Parameters.SlotsToCoeffsParameters.Scaling)
	require.NotNil(t, fastEval.Parameters.SlotsToCoeffsParameters.Scaling)
	require.Zero(t, standardEval.Parameters.SlotsToCoeffsParameters.Scaling.Cmp(fastEval.Parameters.SlotsToCoeffsParameters.Scaling))
}

func TestFastBootstrapOrdinaryStageAPI(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	sk := rlwe.NewKeyGenerator(params.BootstrappingParameters).GenSecretKeyNew()
	keys, _, err := params.GenEvaluationKeys(sk)
	require.NoError(t, err)
	eval, err := NewEvaluator(params, keys)
	require.NoError(t, err)

	input := fastBootstrapPlainCiphertext(t, residual, 0, 2)
	packed, ctxtN1, ctxtN2, err := eval.PackAndSwitchN1ToN2([]rlwe.Ciphertext{*input})
	require.NoError(t, err)
	require.Len(t, packed, 1)

	var scaled *rlwe.Ciphertext
	scaled, _, err = eval.ScaleDown(&packed[0])
	require.NoError(t, err)
	packed[0] = *scaled
	var modded *rlwe.Ciphertext
	modded, err = eval.ModUp(&packed[0])
	require.NoError(t, err)
	packed[0] = *modded
	ctReal, ctImag, err := eval.CoeffsToSlots(&packed[0])
	require.NoError(t, err)
	require.NotNil(t, ctReal)
	ctReal, err = eval.EvalMod(ctReal)
	require.NoError(t, err)
	if ctImag != nil {
		ctImag, err = eval.EvalMod(ctImag)
		require.NoError(t, err)
	}
	ctOut, err := eval.SlotsToCoeffs(ctReal, ctImag)
	require.NoError(t, err)
	outputs, err := eval.UnpackAndSwitchN2ToN1([]rlwe.Ciphertext{*ctOut}, ctxtN1, ctxtN2)
	require.NoError(t, err)
	require.Len(t, outputs, 1)
	require.Equal(t, residual.MaxLevel(), outputs[0].Level())
	require.True(t, outputs[0].IsNTT)
	require.False(t, outputs[0].IsMontgomery)
	require.True(t, outputs[0].Scale.Equal(residual.DefaultScale()))

	full, err := eval.Bootstrap(fastBootstrapPlainCiphertext(t, residual, 0, 2))
	require.NoError(t, err)
	require.Equal(t, outputs[0].Level(), full.Level())
	require.True(t, outputs[0].Scale.Equal(full.Scale))
}

func TestFastBootstrapMatchesStandardDecodedReference(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	fastEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	values := []complex128{0.125 + 0.25i, -0.25 + 0.0625i, 0.375 - 0.125i, -0.0625 - 0.1875i}
	fastInput := fastBootstrapEncodedCiphertext(t, residual, 2, values)
	fastOut, err := fastEval.Bootstrap(fastInput)
	require.NoError(t, err)

	sk := rlwe.NewKeyGenerator(params.BootstrappingParameters).GenSecretKeyNew()
	keys, _, err := params.GenEvaluationKeys(sk)
	require.NoError(t, err)
	standardEval, err := NewEvaluator(params, keys)
	require.NoError(t, err)
	standardInput := fastBootstrapEncodedCiphertext(t, residual, 2, values)
	standardOut, err := standardEval.Bootstrap(standardInput)
	require.NoError(t, err)
	decryptor := rlwe.NewDecryptor(residual, sk)
	standardPlaintext := decryptor.DecryptNew(standardOut)
	standardValues := make([]complex128, len(values))
	require.NoError(t, ckks.NewEncoder(residual).Decode(standardPlaintext, standardValues))
	fastPlaintext := ckks.NewPlaintext(residual, fastOut.Level())
	*fastPlaintext.MetaData = *fastOut.MetaData
	fastPlaintext.Value.Copy(fastOut.Value[0])
	fastPlaintext.IsNTT = fastOut.IsNTT
	fastPlaintext.IsMontgomery = fastOut.IsMontgomery
	fastValues := make([]complex128, len(values))
	require.NoError(t, ckks.NewEncoder(residual).Decode(fastPlaintext, fastValues))
	for i := range values {
		require.InDelta(t, real(standardValues[i]), real(fastValues[i]), 1e-2, "real slot %d", i)
		require.InDelta(t, imag(standardValues[i]), imag(fastValues[i]), 1e-2, "imag slot %d", i)
	}
}

func TestFastBootstrapCoreIgnoresDormantResidues(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.GreaterOrEqual(t, params.BootstrappingParameters.MaxLevel(), 2)
	first := fastBootstrapPlainCiphertext(t, params.BootstrappingParameters, 2, 2)
	q01Scale := new(big.Float).SetUint64(params.BootstrappingParameters.Q()[0])
	q01Scale.Mul(q01Scale, new(big.Float).SetUint64(params.BootstrappingParameters.Q()[1]))
	first.Scale = rlwe.NewScale(q01Scale)
	second := first.CopyNew()
	for d := range second.Value {
		for limb := 2; limb < len(second.Value[d].Coeffs); limb++ {
			for i := range second.Value[d].Coeffs[limb] {
				second.Value[d].Coeffs[limb][i] = uint64(i + 17 + d + limb)
			}
		}
	}
	firstOut, _, err := eval.bootstrapCore(first, 0)
	require.NoError(t, err)
	secondOut, _, err := eval.bootstrapCore(second, 0)
	require.NoError(t, err)
	for d := 0; d <= 1; d++ {
		for limb := 0; limb <= secondOut.Level() && limb < 2; limb++ {
			require.Equal(t, firstOut.Value[d].Coeffs[limb], secondOut.Value[d].Coeffs[limb], "component %d limb %d", d, limb)
		}
	}
}

func TestFastBootstrapCoreS2CConsumesEvalModQ3(t *testing.T) {
	params, residual := fastBootstrapParameters(t, 2, false)
	params.ResidualParameters = residual
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	require.NoError(t, eval.ensureFastBootstrapCircuit())
	require.Equal(t, 3, eval.S2CDFTMatrix.LevelQ)
	require.Equal(t, eval.S2CDFTMatrix.LevelQ, eval.Mod1Parameters.LevelQ-eval.Parameters.Mod1ParametersLiteral.Depth())

	input := fastBootstrapPlainCiphertext(t, residual, residual.MaxLevel(), 2)
	baseline, _, err := eval.bootstrapCore(input.CopyNew(), 0)
	require.NoError(t, err)
	require.Equal(t, 1, baseline.Level(), "two S2C groups should contract Level 3 to Level 1")

	// Poison only q3 in the first production S2C factor. The real bootstrapCore
	// path must have passed all four EvalMod rows to S2C, so this change must
	// affect the q012 result after the first Rescale contracts Level 3 to 2.
	changedQ3 := 0
	for diagonal, poly := range eval.S2CDFTMatrix.Matrices[0].Vec {
		if len(poly.Q.Coeffs) <= 3 {
			continue
		}
		for i, value := range poly.Q.Coeffs[3] {
			if value != 0 {
				changedQ3++
			}
			poly.Q.Coeffs[3][i] = 0
		}
		eval.S2CDFTMatrix.Matrices[0].Vec[diagonal] = poly
	}
	require.Greater(t, changedQ3, 0, "first S2C factor must contain encoded q3 data")

	poisoned, _, err := eval.bootstrapCore(input.CopyNew(), 0)
	require.NoError(t, err)
	require.Equal(t, baseline.Level(), poisoned.Level())
	require.Equal(t, baseline.Scale, poisoned.Scale)
	different := false
	for component := range baseline.Value {
		for row := 0; row <= baseline.Level(); row++ {
			for i, value := range baseline.Value[component].Coeffs[row] {
				if value != poisoned.Value[component].Coeffs[row][i] {
					different = true
					break
				}
			}
			if different {
				break
			}
		}
		if different {
			break
		}
	}
	require.True(t, different, "production Bootstrap S2C output must depend on the authoritative EvalMod q3 row")
}
