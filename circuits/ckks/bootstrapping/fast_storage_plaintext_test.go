package bootstrapping

import (
	"math/big"
	"runtime"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	ckkslintrans "github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	commonlintrans "github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

func TestFASTSTORAGE007ActualLogN13C2SFeasibility(t *testing.T) {
	params, _ := fastLogN13CompressionParameters(t)
	bootstrapParams := params.BootstrappingParameters
	_, _, logicalC2S, _, err := buildBootstrapCircuitData(params)
	require.NoError(t, err)
	fastC2S, restorePlan, active, err := prepareFastLogN13C2S(params, logicalC2S)
	require.NoError(t, err)
	require.True(t, active)
	require.Len(t, restorePlan, len(fastC2S.Levels))

	// Reuse an existing deterministic bounded Fast-Bootstrap input fixture,
	// import its exact level-0 centered lifts, and apply the existing private-F
	// ModUp boundary. This yields reproducible per-component input bounds.
	source := fastBootstrapPlainCiphertext(t, bootstrapParams, 0, params.LogMaxSlots())
	private, err := fastckks.ImportLevel0(bootstrapParams, source, 3, fastckks.FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	private, err = fastckks.FastStorageModUpLevel0(private, bootstrapParams.MaxLevel())
	require.NoError(t, err)
	currentBounds := private.ComponentBounds()
	require.Len(t, currentBounds, 2)

	s3, err := fastckks.FastStorageProduct(bootstrapParams, 3)
	require.NoError(t, err)

	var representative *fastckks.FastStorageLinearTransformation
	var representativeLogical ckkslintrans.LinearTransformation
	for factor, logicalMatrix := range fastC2S.Matrices {
		runtime.GC()
		var before, after runtime.MemStats
		runtime.ReadMemStats(&before)
		started := time.Now()
		mirrored, err := fastckks.NewFastStorageLinearTransformation(bootstrapParams, commonlintrans.LinearTransformation(logicalMatrix))
		require.NoError(t, err, "mirror C2S factor %d", factor)
		buildTime := time.Since(started)
		runtime.ReadMemStats(&after)
		retainedBytes := uint64(mirrored.DiagonalCount() * 3 * bootstrapParams.N() * 8)
		allocatedBytes := after.TotalAlloc - before.TotalAlloc
		matrixScale := mirrored.Scale()
		lSigma := mirrored.SumDiagonalL1Norm()
		factorBounds := make([]*big.Int, len(currentBounds))
		fits := true
		for component, bound := range currentBounds {
			factorBounds[component] = new(big.Int).Mul(bound, lSigma)
			doubled := new(big.Int).Lsh(new(big.Int).Set(factorBounds[component]), 1)
			if doubled.Cmp(s3) >= 0 {
				fits = false
			}
		}
		t.Logf("FAST-STORAGE-007 factor=%d levelQ=%d scale=%s diagonals=%d N1=%d bsgs=%t max_Bp=%s max_Lp=%s Lsigma=%s input_bounds=[%s,%s] output_bounds=[%s,%s] strict_width3=%t build=%s retained_bytes=%d total_alloc_bytes=%d",
			factor, mirrored.LevelQ(), matrixScale.Value.Text('g', 12), mirrored.DiagonalCount(), mirrored.N1(), mirrored.N1() != 0,
			mirrored.MaxAbsCoefficient(), mirrored.MaxDiagonalL1Norm(), lSigma,
			currentBounds[0], currentBounds[1], factorBounds[0], factorBounds[1], fits,
			buildTime, retainedBytes, allocatedBytes)
		if factor == 0 {
			representative = mirrored
			representativeLogical = logicalMatrix
		}
		if !fits {
			t.Logf("classification=FAST_STORAGE_007_FOUNDATION_VALID_C2S_CAPACITY_BLOCKED first_failing_factor=%d S3=%s doubled_output_bounds=[%s,%s]", factor, s3, new(big.Int).Lsh(new(big.Int).Set(factorBounds[0]), 1), new(big.Int).Lsh(new(big.Int).Set(factorBounds[1]), 1))
			return
		}

		// Propagate the conceptual C2S group bound with the frozen 004 logical
		// Rescale theorem, then apply the already accepted restore scalar.
		q := bootstrapParams.Q()[private.LogicalLevel()]
		for component, bound := range factorBounds {
			numerator := new(big.Int).Add(bound, new(big.Int).SetUint64((q-1)/2))
			currentBounds[component] = numerator.Quo(numerator, new(big.Int).SetUint64(q))
		}
		if err := private.SetLogicalLevel(private.LogicalLevel() - 1); err != nil {
			require.NoError(t, err)
		}
		if restorePlan[factor] != 0 {
			scalar := new(big.Int).Lsh(big.NewInt(1), uint(restorePlan[factor]))
			for component := range currentBounds {
				currentBounds[component].Mul(currentBounds[component], scalar)
			}
		}
		t.Logf("FAST-STORAGE-007 group=%d rescale_q_index=%d rescale_q=%d restore_exp=%d propagated_bounds=[%s,%s]", factor, private.LogicalLevel()+1, q, restorePlan[factor], currentBounds[0], currentBounds[1])
	}
	t.Log("classification=FAST_STORAGE_007_PRIVATE_F_LINEAR_TRANSFORM_FEASIBLE_CANDIDATE")
	require.NotNil(t, representative)

	// Compare both arithmetic domains from one identical level-16 bounded
	// private-F input. Convert q rows only at the explicit export boundary;
	// Montgomery residues are never reinterpreted as ordinary residues.
	privateInput, err := fastckks.ImportLevel0(bootstrapParams, source, 3, fastckks.FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	privateInput, err = fastckks.FastStorageModUpLevel0(privateInput, bootstrapParams.MaxLevel())
	require.NoError(t, err)
	privateOutput, err := fastckks.FastStorageLinearTransform(privateInput, representative)
	require.NoError(t, err)
	privateLogical, err := privateOutput.ExportToCompactLogical(fastckks.FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)

	logicalInput, err := privateInput.ExportToCompactLogical(fastckks.FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	logicalRing := bootstrapParams.RingQ().AtLevel(logicalInput.Level())
	maintained := fastckks.MaintainedLimbCount(&bootstrapParams, logicalInput.Level())
	for component := 0; component < 2; component++ {
		for row := 0; row < maintained; row++ {
			logicalRing.SubRings[row].MForm(logicalInput.Value[component].Coeffs[row], logicalInput.Value[component].Coeffs[row])
		}
	}
	logicalInput.IsMontgomery = true
	logicalOutput := ckks.NewCiphertext(bootstrapParams, 1, logicalInput.Level())
	logicalOutput.IsNTT, logicalOutput.IsMontgomery = true, true
	fastEval := fastckks.NewEvaluator(bootstrapParams)
	require.NoError(t, fastEval.LinearTransform(logicalInput, commonlintrans.LinearTransformation(representativeLogical), logicalOutput))
	require.Equal(t, privateLogical.Level(), logicalOutput.Level())
	require.True(t, privateLogical.Scale.Equal(logicalOutput.Scale))
	for component := 0; component < 2; component++ {
		for row := 0; row < maintained; row++ {
			logicalRing.SubRings[row].IMForm(logicalOutput.Value[component].Coeffs[row], logicalOutput.Value[component].Coeffs[row])
		}
	}
	for component := range privateLogical.Value {
		for row := 0; row < maintained; row++ {
			require.Equal(t, privateLogical.Value[component].Coeffs[row], logicalOutput.Value[component].Coeffs[row], "factor 0 component=%d q%d", component, row)
		}
	}
	outputScale := privateOutput.Scale()
	t.Logf("FAST-STORAGE-007 factor0_exact_logical_fast_match=true maintained_q_rows=%d output_scale=%s", maintained, outputScale.Value.Text('g', 12))
}

func BenchmarkFASTSTORAGE007LogN13PrivateFLinearTransform(b *testing.B) {
	params, _ := fastLogN13CompressionParameters(b)
	bootstrapParams := params.BootstrappingParameters
	_, _, logicalC2S, _, err := buildBootstrapCircuitData(params)
	if err != nil {
		b.Fatal(err)
	}
	fastC2S, _, active, err := prepareFastLogN13C2S(params, logicalC2S)
	if err != nil {
		b.Fatal(err)
	}
	if !active {
		b.Fatal("expected actual LogN13 Fast C2S preparation")
	}
	matrix := fastC2S.Matrices[0]
	source := fastBootstrapPlainCiphertext(b, bootstrapParams, 0, params.LogMaxSlots())
	private, err := fastckks.ImportLevel0(bootstrapParams, source, 3, fastckks.FastCiphertextDomain{IsNTT: true})
	if err != nil {
		b.Fatal(err)
	}
	private, err = fastckks.FastStorageModUpLevel0(private, bootstrapParams.MaxLevel())
	if err != nil {
		b.Fatal(err)
	}
	privateMatrix, err := fastckks.NewFastStorageLinearTransformation(bootstrapParams, commonlintrans.LinearTransformation(matrix))
	if err != nil {
		b.Fatal(err)
	}
	logical, err := private.ExportToCompactLogical(fastckks.FastCiphertextDomain{IsNTT: true})
	if err != nil {
		b.Fatal(err)
	}
	logicalRing := bootstrapParams.RingQ().AtLevel(logical.Level())
	maintained := fastckks.MaintainedLimbCount(&bootstrapParams, logical.Level())
	for component := 0; component < 2; component++ {
		for row := 0; row < maintained; row++ {
			logicalRing.SubRings[row].MForm(logical.Value[component].Coeffs[row], logical.Value[component].Coeffs[row])
		}
	}
	logical.IsMontgomery = true
	logicalOut := ckks.NewCiphertext(bootstrapParams, 1, logical.Level())
	logicalOut.IsNTT, logicalOut.IsMontgomery = true, true
	eval := fastckks.NewEvaluator(bootstrapParams)
	b.ReportMetric(float64(privateMatrix.DiagonalCount()), "diagonals")
	b.Run("PrivateFWidth3", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, err := fastckks.FastStorageLinearTransform(private, privateMatrix); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("LogicalFast", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := eval.LinearTransform(logical, commonlintrans.LinearTransformation(matrix), logicalOut); err != nil {
				b.Fatal(err)
			}
		}
	})
}
