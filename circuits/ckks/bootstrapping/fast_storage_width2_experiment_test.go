package bootstrapping

import (
	"fmt"
	"math/big"
	"runtime"
	"strings"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	commonlintrans "github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

type fastStorage008Fixture struct {
	params        Parameters
	bootstrap     ckks.Parameters
	matrices      []commonlintrans.LinearTransformation
	restorePlan   []int
	privateWidth3 *fastckks.FastCiphertext
	privateWidth2 *fastckks.FastCiphertext
}

func newFastStorage008Fixture(tb testing.TB) fastStorage008Fixture {
	tb.Helper()
	params, _ := fastLogN13CompressionParameters(tb)
	bootstrap := params.BootstrappingParameters
	_, _, logicalC2S, _, err := buildBootstrapCircuitData(params)
	require.NoError(tb, err)
	fastC2S, restorePlan, active, err := prepareFastLogN13C2S(params, logicalC2S)
	require.NoError(tb, err)
	require.True(tb, active, "FAST-STORAGE-008 requires the actual compressed LogN13 C2S profile")
	require.Len(tb, fastC2S.Matrices, 4)
	require.Equal(tb, []int{4, 2, 0, 0}, restorePlan)

	matrices := make([]commonlintrans.LinearTransformation, len(fastC2S.Matrices))
	for i := range fastC2S.Matrices {
		matrices[i] = commonlintrans.LinearTransformation(fastC2S.Matrices[i])
	}
	source := fastBootstrapPlainCiphertext(tb, bootstrap, 0, params.LogMaxSlots())
	private3, err := fastckks.ImportLevel0(bootstrap, source, 3, fastckks.FastCiphertextDomain{IsNTT: true})
	require.NoError(tb, err)
	private3, err = fastckks.FastStorageModUpLevel0(private3, bootstrap.MaxLevel())
	require.NoError(tb, err)
	private2, err := private3.ContractStorage(2)
	require.NoError(tb, err, "same ModUp input must fit width 2 before the C2S chain")
	return fastStorage008Fixture{
		params: params, bootstrap: bootstrap, matrices: matrices,
		restorePlan: restorePlan, privateWidth3: private3, privateWidth2: private2,
	}
}

func TestFASTSTORAGE008ActualLogN13C2SWidth2Experiment(t *testing.T) {
	fixture := newFastStorage008Fixture(t)
	s2, err := fastckks.FastStorageProduct(fixture.bootstrap, 2)
	require.NoError(t, err)
	s3, err := fastckks.FastStorageProduct(fixture.bootstrap, 3)
	require.NoError(t, err)
	t.Logf("FAST-STORAGE-008 LogN=%d N=%d S2=%s S3=%s input_bounds=%v", fixture.bootstrap.LogN(), fixture.bootstrap.N(), s2, s3, fixture.privateWidth3.ComponentBounds())

	current3 := fixture.privateWidth3.CopyNew()
	current2 := fixture.privateWidth2.CopyNew()
	logical := fastStorage008LogicalInput(t, fixture.bootstrap, current3)
	logicalEval := fastckks.NewEvaluator(fixture.bootstrap)
	width2Blocked := false
	classification := "FAST_STORAGE_008_WIDTH2_VALID_NOT_COMPETITIVE"

	for factor, logicalMatrix := range fixture.matrices {
		var transformBounds2, groupBounds2 []*big.Int
		mirror3, build3, alloc3 := fastStorage008BuildMirror(t, fixture.bootstrap, logicalMatrix, 3)
		mirror2, build2, alloc2, err2 := fastStorage008BuildMirrorMaybe(t, fixture.bootstrap, logicalMatrix, 2)
		if err2 != nil {
			require.True(t, isFastStorage008CapacityFailure(err2), "unexpected width-2 mirror error: %v", err2)
			width2Blocked = true
			classification = "FAST_STORAGE_008_WIDTH2_CAPACITY_BLOCKED"
			t.Logf("factor=%d stage=plaintext_mirror width=2 classification=%s error=%q", factor, classification, err2)
		}
		boundsIn3 := current3.ComponentBounds()
		fastStorage008RequireCapacity(t, boundsIn3, s3, 3, factor, "input")
		out3, err := fastckks.FastStorageLinearTransform(current3, mirror3)
		require.NoError(t, err, "width-3 transform factor %d must fit", factor)
		fastStorage008RequireCapacity(t, out3.ComponentBounds(), s3, 3, factor, "linear_transform")

		logicalOut := ckks.NewCiphertext(fixture.bootstrap, 1, logical.Level())
		require.NoError(t, logicalEval.LinearTransform(logical, logicalMatrix, logicalOut), "Logical Fast factor %d", factor)
		require.Equal(t, out3.LogicalLevel(), logicalOut.Level())
		require.True(t, out3.Scale().Equal(logicalOut.Scale), "factor %d transform Scale", factor)

		rescaled3, err := fastckks.FastStorageRescale(out3)
		require.NoError(t, err, "width-3 Rescale factor %d", factor)
		fastStorage008RequireCapacity(t, rescaled3.ComponentBounds(), s3, 3, factor, "rescale")
		if exponent := fixture.restorePlan[factor]; exponent != 0 {
			rescaled3, err = fastckks.FastStorageMulInteger(rescaled3, new(big.Int).Lsh(big.NewInt(1), uint(exponent)))
			require.NoError(t, err, "width-3 restore factor %d", factor)
			fastStorage008RequireCapacity(t, rescaled3.ComponentBounds(), s3, 3, factor, "restore")
		}

		logicalRescaled := ckks.NewCiphertext(fixture.bootstrap, 1, logicalOut.Level()-1)
		require.NoError(t, logicalEval.Rescale(logicalOut, logicalRescaled), "Logical Fast Rescale factor %d", factor)
		if exponent := fixture.restorePlan[factor]; exponent != 0 {
			logicalRestored := ckks.NewCiphertext(fixture.bootstrap, 1, logicalRescaled.Level())
			logicalRestored.IsNTT = logicalRescaled.IsNTT
			logicalRestored.IsMontgomery = logicalRescaled.IsMontgomery
			require.NoError(t, logicalEval.MulIntegerMaintained(logicalRescaled, new(big.Int).Lsh(big.NewInt(1), uint(exponent)), logicalRestored))
			logicalRescaled = logicalRestored
		}
		fastStorage008LogicalRowsEqual(t, fixture.bootstrap, rescaled3, logicalRescaled, fmt.Sprintf("factor %d width-3", factor))

		if !width2Blocked {
			require.NotNil(t, mirror2)
			fastStorage008RequireCapacity(t, current2.ComponentBounds(), s2, 2, factor, "input")
			out2, transformErr := fastckks.FastStorageLinearTransform(current2, mirror2)
			if transformErr != nil {
				require.True(t, isFastStorage008CapacityFailure(transformErr), "unexpected width-2 transform error at factor %d: %v", factor, transformErr)
				width2Blocked = true
				classification = "FAST_STORAGE_008_WIDTH2_CAPACITY_BLOCKED"
				t.Logf("factor=%d stage=linear_transform width=2 classification=%s input_bounds=%v bound=%v error=%q", factor, classification, current2.ComponentBounds(), mirror2.SumDiagonalL1Norm(), transformErr)
			} else {
				transformBounds2 = out2.ComponentBounds()
				fastStorage008RequireCapacity(t, out2.ComponentBounds(), s2, 2, factor, "linear_transform")
				rescaled2, rescaleErr := fastckks.FastStorageRescale(out2)
				require.NoError(t, rescaleErr, "width-2 Rescale factor %d", factor)
				fastStorage008RequireCapacity(t, rescaled2.ComponentBounds(), s2, 2, factor, "rescale")
				if exponent := fixture.restorePlan[factor]; exponent != 0 {
					rescaled2, rescaleErr = fastckks.FastStorageMulInteger(rescaled2, new(big.Int).Lsh(big.NewInt(1), uint(exponent)))
					require.NoError(t, rescaleErr, "width-2 restore factor %d", factor)
					fastStorage008RequireCapacity(t, rescaled2.ComponentBounds(), s2, 2, factor, "restore")
				}
				equalLift, compareErr := fastckks.FastStorageEqualLift(rescaled2, rescaled3)
				require.NoError(t, compareErr)
				require.True(t, equalLift, "width-2 and width-3 centered lifts differ after factor %d", factor)
				fastStorage008LogicalRowsEqual(t, fixture.bootstrap, rescaled2, logicalRescaled, fmt.Sprintf("factor %d width-2", factor))
				current2 = rescaled2
				groupBounds2 = rescaled2.ComponentBounds()
			}
		}

		t.Logf("factor=%d diagonals=%d N1=%d build_width2=%s/%dB build_width3=%s/%dB input_bounds=%v transform_bounds_width3=%v group_bounds_width3=%v transform_bounds_width2=%v group_bounds_width2=%v width2_capacity=%t exact_width2_width3_lift=%t logical_rows_width3=true",
			factor, mirror3.DiagonalCount(), mirror3.N1(), build2, alloc2, build3, alloc3,
			boundsIn3, out3.ComponentBounds(), rescaled3.ComponentBounds(), transformBounds2, groupBounds2, !width2Blocked, !width2Blocked)
		current3 = rescaled3
		logical = logicalRescaled
	}

	if width2Blocked {
		t.Logf("classification=%s; width-2 chain stopped at first capacity failure; width-3 and Logical Fast completed all four factors", classification)
		return
	}
	require.NotNil(t, current2)
	require.Equal(t, "FAST_STORAGE_008_WIDTH2_VALID_NOT_COMPETITIVE", classification)
	t.Logf("classification=%s; all four width-2 stages passed strict capacity and exact semantic checks", classification)
}

func fastStorage008BuildMirror(t testing.TB, params ckks.Parameters, matrix commonlintrans.LinearTransformation, width int) (*fastckks.FastStorageLinearTransformation, time.Duration, uint64) {
	t.Helper()
	mirrored, build, allocated, err := fastStorage008BuildMirrorMaybe(t, params, matrix, width)
	require.NoError(t, err, "mirror width-%d C2S matrix", width)
	return mirrored, build, allocated
}

func fastStorage008BuildMirrorMaybe(t testing.TB, params ckks.Parameters, matrix commonlintrans.LinearTransformation, width int) (*fastckks.FastStorageLinearTransformation, time.Duration, uint64, error) {
	t.Helper()
	runtime.GC()
	var before, after runtime.MemStats
	runtime.ReadMemStats(&before)
	started := time.Now()
	mirrored, err := fastckks.NewFastStorageLinearTransformationWithWidth(params, matrix, width)
	build := time.Since(started)
	runtime.ReadMemStats(&after)
	return mirrored, build, after.TotalAlloc - before.TotalAlloc, err
}

func isFastStorage008CapacityFailure(err error) bool {
	message := strings.ToLower(err.Error())
	return strings.Contains(message, "capacity") || strings.Contains(message, "centered")
}

func fastStorage008RequireCapacity(t *testing.T, bounds []*big.Int, storageProduct *big.Int, width, factor int, stage string) {
	t.Helper()
	for component, bound := range bounds {
		doubled := new(big.Int).Lsh(new(big.Int).Set(bound), 1)
		require.Less(t, doubled.Cmp(storageProduct), 0, "width-%d factor=%d stage=%s component=%d strict 2B<S%d failed: bound=%s S=%s", width, factor, stage, component, width, bound, storageProduct)
	}
}

func fastStorage008LogicalInput(t testing.TB, params ckks.Parameters, private *fastckks.FastCiphertext) *rlwe.Ciphertext {
	t.Helper()
	logical, err := private.ExportToCompactLogical(fastckks.FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	fastStorage008SetMontgomery(params, logical)
	return logical
}

func fastStorage008SetMontgomery(params ckks.Parameters, ct *rlwe.Ciphertext) {
	maintained := fastckks.MaintainedLimbCount(&params, ct.Level())
	ringQ := params.RingQ().AtLevel(ct.Level())
	for component := range ct.Value {
		for row := 0; row < maintained; row++ {
			ringQ.SubRings[row].MForm(ct.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
		}
	}
	ct.IsMontgomery = true
}

func fastStorage008LogicalRowsEqual(t *testing.T, params ckks.Parameters, private *fastckks.FastCiphertext, logicalMontgomery *rlwe.Ciphertext, context string) {
	t.Helper()
	got, err := private.ExportToCompactLogical(fastckks.FastCiphertextDomain{IsNTT: true})
	require.NoError(t, err)
	want := logicalMontgomery.CopyNew()
	maintained := fastckks.MaintainedLimbCount(&params, got.Level())
	ringQ := params.RingQ().AtLevel(got.Level())
	for component := range want.Value {
		for row := 0; row < maintained; row++ {
			ringQ.SubRings[row].IMForm(want.Value[component].Coeffs[row], want.Value[component].Coeffs[row])
		}
	}
	want.IsMontgomery = false
	require.Equal(t, got.Level(), want.Level(), "%s logical level", context)
	require.True(t, got.Scale.Equal(want.Scale), "%s logical scale", context)
	for component := range got.Value {
		for row := 0; row < maintained; row++ {
			require.Equal(t, got.Value[component].Coeffs[row], want.Value[component].Coeffs[row], "%s component=%d q%d", context, component, row)
		}
	}
}

func fastStorage008PrivateChain(input *fastckks.FastCiphertext, matrices []*fastckks.FastStorageLinearTransformation, restore []int) (*fastckks.FastCiphertext, error) {
	current := input.CopyNew()
	for i, matrix := range matrices {
		var err error
		current, err = fastckks.FastStorageLinearTransform(current, matrix)
		if err != nil {
			return nil, fmt.Errorf("factor %d LinearTransform: %w", i, err)
		}
		current, err = fastckks.FastStorageRescale(current)
		if err != nil {
			return nil, fmt.Errorf("factor %d Rescale: %w", i, err)
		}
		if restore[i] != 0 {
			current, err = fastckks.FastStorageMulInteger(current, new(big.Int).Lsh(big.NewInt(1), uint(restore[i])))
			if err != nil {
				return nil, fmt.Errorf("factor %d restore: %w", i, err)
			}
		}
	}
	return current, nil
}

func fastStorage008LogicalChain(input *rlwe.Ciphertext, params ckks.Parameters, eval *fastckks.Evaluator, matrices []commonlintrans.LinearTransformation, restore []int) (*rlwe.Ciphertext, error) {
	current := input.CopyNew()
	for i, matrix := range matrices {
		transformed := ckks.NewCiphertext(params, 1, current.Level())
		if err := eval.LinearTransform(current, matrix, transformed); err != nil {
			return nil, fmt.Errorf("factor %d LinearTransform: %w", i, err)
		}
		rescaled := ckks.NewCiphertext(params, 1, transformed.Level()-1)
		if err := eval.Rescale(transformed, rescaled); err != nil {
			return nil, fmt.Errorf("factor %d Rescale: %w", i, err)
		}
		if restore[i] != 0 {
			restored := ckks.NewCiphertext(params, 1, rescaled.Level())
			restored.IsNTT = rescaled.IsNTT
			restored.IsMontgomery = rescaled.IsMontgomery
			if err := eval.MulIntegerMaintained(rescaled, new(big.Int).Lsh(big.NewInt(1), uint(restore[i])), restored); err != nil {
				return nil, fmt.Errorf("factor %d restore: %w", i, err)
			}
			rescaled = restored
		}
		current = rescaled
	}
	return current, nil
}

func BenchmarkFASTSTORAGE008LogN13C2SFactor0(b *testing.B) {
	fixture := newFastStorage008Fixture(b)
	matrix3, _, _ := fastStorage008BuildMirror(b, fixture.bootstrap, fixture.matrices[0], 3)
	matrix2, _, _, err2 := fastStorage008BuildMirrorMaybe(b, fixture.bootstrap, fixture.matrices[0], 2)
	logicalInput := fastStorage008LogicalInput(b, fixture.bootstrap, fixture.privateWidth3)
	logicalOut := ckks.NewCiphertext(fixture.bootstrap, 1, logicalInput.Level())
	logicalEval := fastckks.NewEvaluator(fixture.bootstrap)
	b.Run("LogicalFast", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if err := logicalEval.LinearTransform(logicalInput, fixture.matrices[0], logicalOut); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("PrivateFWidth3", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, err := fastckks.FastStorageLinearTransform(fixture.privateWidth3, matrix3); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("PrivateFWidth2", func(b *testing.B) {
		if err2 != nil {
			b.Skipf("width-2 factor-0 mirror capacity blocked: %v", err2)
		}
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, err := fastckks.FastStorageLinearTransform(fixture.privateWidth2, matrix2); err != nil {
				b.Fatal(err)
			}
		}
	})
}

func BenchmarkFASTSTORAGE008LogN13C2SFullChain(b *testing.B) {
	fixture := newFastStorage008Fixture(b)
	matrices3 := make([]*fastckks.FastStorageLinearTransformation, len(fixture.matrices))
	matrices2 := make([]*fastckks.FastStorageLinearTransformation, len(fixture.matrices))
	var width2Err error
	for i, matrix := range fixture.matrices {
		matrices3[i], _, _ = fastStorage008BuildMirror(b, fixture.bootstrap, matrix, 3)
		if width2Err == nil {
			matrices2[i], _, _, width2Err = fastStorage008BuildMirrorMaybe(b, fixture.bootstrap, matrix, 2)
		}
	}
	logicalInput := fastStorage008LogicalInput(b, fixture.bootstrap, fixture.privateWidth3)
	logicalEval := fastckks.NewEvaluator(fixture.bootstrap)
	b.Run("LogicalFast", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, err := fastStorage008LogicalChain(logicalInput, fixture.bootstrap, logicalEval, fixture.matrices, fixture.restorePlan); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("PrivateFWidth3", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, err := fastStorage008PrivateChain(fixture.privateWidth3, matrices3, fixture.restorePlan); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("PrivateFWidth2", func(b *testing.B) {
		if width2Err != nil {
			b.Skipf("width-2 chain capacity blocked: %v", width2Err)
		}
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, err := fastStorage008PrivateChain(fixture.privateWidth2, matrices2, fixture.restorePlan); err != nil {
				b.Fatal(err)
			}
		}
	})
}
