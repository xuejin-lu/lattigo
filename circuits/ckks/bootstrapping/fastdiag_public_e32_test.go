//go:build fastdiag

package bootstrapping

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"strconv"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/internal/fastdiag"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

type fastDiagPublicE32Manifest struct {
	SchemaVersion               string              `json:"schema_version"`
	Profile                     string              `json:"profile"`
	Backend                     string              `json:"backend"`
	BackendCommit               string              `json:"backend_commit"`
	PrimaryCommit               string              `json:"primary_commit"`
	PrimarySourceSHA256         string              `json:"primary_measurement_source_sha256"`
	ConfigSHA256                string              `json:"config_sha256"`
	QPSHA256                    string              `json:"qp_sha256"`
	InputSHA256                 string              `json:"canonical_input_sha256"`
	WorkloadSHA256              string              `json:"workload_sha256"`
	ParametersFile              string              `json:"parameters_file"`
	ParametersSHA256            string              `json:"parameters_sha256"`
	CiphertextFile              string              `json:"ciphertext_file"`
	CiphertextSHA256            string              `json:"ciphertext_sha256"`
	CiphertextFingerprintSHA256 string              `json:"ciphertext_fingerprint_sha256"`
	DecodedSHA256               string              `json:"decoded_sha256"`
	State                       fastDiagPublicState `json:"ciphertext_state"`
	PhysicalRowLengths          [][]int             `json:"physical_row_lengths"`
	IsNTT                       bool                `json:"is_ntt"`
	IsMontgomery                bool                `json:"is_montgomery"`
	C1Nonzero                   bool                `json:"c1_nonzero"`
	NumericalGate               float64             `json:"max_complex_numerical_gate"`
	EphemeralSecretWeight       int                 `json:"ephemeral_secret_weight"`
	Parameters                  fastDiagEffective   `json:"effective_parameters"`
}

type fastDiagPublicState struct {
	Level       int    `json:"level"`
	Degree      int    `json:"degree"`
	Scale       string `json:"scale"`
	QPrefixRows int    `json:"q_prefix_rows"`
	PrefixQ     string `json:"prefix_q"`
}

type fastDiagEffective struct {
	LogN         int      `json:"log_n"`
	LogSlots     int      `json:"log_slots"`
	InputSlots   int      `json:"input_slots"`
	RingN        int      `json:"ring_n"`
	QPrimes      []string `json:"q_primes"`
	PPrimes      []string `json:"p_primes"`
	DefaultScale string   `json:"default_scale"`
	Mod1Scale    int      `json:"mod1_log_scale"`
	Mod1Degree   int      `json:"mod1_degree"`
	DoubleAngle  int      `json:"double_angle"`
	K            int      `json:"k"`
	CircuitOrder int      `json:"circuit_order"`
}

type fastDiagPublicVector struct {
	Real float64 `json:"real"`
	Imag float64 `json:"imag"`
}

type fastDiagPublicVectors struct {
	SchemaVersion    string                            `json:"schema_version"`
	Mode             string                            `json:"mode"`
	Backend          string                            `json:"backend"`
	BackendCommit    string                            `json:"backend_commit"`
	PrimaryCommit    string                            `json:"primary_commit"`
	SourceSHA256     string                            `json:"primary_measurement_source_sha256"`
	ConfigSHA256     string                            `json:"config_sha256"`
	QPSHA256         string                            `json:"qp_sha256"`
	InputSHA256      string                            `json:"input_sha256"`
	WorkloadSHA256   string                            `json:"workload_sha256"`
	Checkpoints      map[string][]fastDiagPublicVector `json:"checkpoints"`
	BootstrapOutputs map[string][]fastDiagPublicVector `json:"bootstrap_outputs"`
}

type fastDiagPublicE32Run struct {
	Index                  int     `json:"index"`
	Phase                  string  `json:"phase"`
	ElapsedNS              int64   `json:"elapsed_ns"`
	DecodedSHA256          string  `json:"decoded_sha256"`
	OracleRMSE             float64 `json:"oracle_complex_rmse"`
	OracleMaxComplex       float64 `json:"oracle_max_complex_difference"`
	FastReferenceWorstRMSE float64 `json:"fast_reference_worst_rmse"`
	FastReferenceWorstMax  float64 `json:"fast_reference_worst_max_complex_difference"`
	Level                  int     `json:"level"`
	Degree                 int     `json:"degree"`
	Scale                  string  `json:"scale"`
	QPrefixRows            int     `json:"q_prefix_rows"`
}

type fastDiagPublicMetrics struct {
	rmse float64
	max  float64
}

type fastDiagPublicE32Artifact struct {
	SchemaVersion         string                 `json:"schema_version"`
	Timestamp             time.Time              `json:"timestamp"`
	TraceStatus           string                 `json:"trace_status"`
	TraceValidationError  string                 `json:"trace_validation_error,omitempty"`
	RawEvidenceFile       string                 `json:"raw_unverified_file,omitempty"`
	RawEvidenceSHA256     string                 `json:"raw_unverified_sha256,omitempty"`
	Profile               string                 `json:"profile"`
	Mode                  string                 `json:"mode"`
	PrimaryCommit         string                 `json:"primary_commit"`
	FixturePrimaryCommit  string                 `json:"fixture_primary_commit"`
	PrimarySourceSHA256   string                 `json:"primary_measurement_source_sha256"`
	ProductionFastCommit  string                 `json:"production_fast_commit"`
	DiagnosticHead        string                 `json:"diagnostic_head"`
	ConfigSHA256          string                 `json:"config_sha256"`
	QPSHA256              string                 `json:"qp_sha256"`
	InputSHA256           string                 `json:"input_sha256"`
	WorkloadSHA256        string                 `json:"workload_sha256"`
	FixtureManifestSHA256 string                 `json:"fixture_manifest_sha256"`
	CiphertextSHA256      string                 `json:"ciphertext_sha256"`
	ExpectedInputSHA256   string                 `json:"expected_input_decoded_sha256"`
	CallBudget            int                    `json:"call_budget"`
	ActualCalls           int                    `json:"actual_calls"`
	CPUProfileFile        string                 `json:"warm_cpu_profile_file"`
	CPUProfileSHA256      string                 `json:"warm_cpu_profile_sha256"`
	HeapProfileFile       string                 `json:"post_warm_heap_profile_file"`
	HeapProfileSHA256     string                 `json:"post_warm_heap_profile_sha256"`
	NumericalGate         float64                `json:"max_complex_gate"`
	GoVersion             string                 `json:"go_version"`
	OS                    string                 `json:"os"`
	Arch                  string                 `json:"arch"`
	NumCPU                int                    `json:"num_cpu"`
	GOMAXPROCS            int                    `json:"gomaxprocs"`
	GOGC                  string                 `json:"gogc"`
	GOMEMLIMIT            string                 `json:"gomemlimit"`
	GODEBUG               string                 `json:"godebug"`
	Runs                  []fastDiagPublicE32Run `json:"runs"`
	Events                []fastdiag.Event       `json:"warm_traced_events"`
}

type fastDiagPublicE32FailureArtifact struct {
	SchemaVersion        string    `json:"schema_version"`
	Timestamp            time.Time `json:"timestamp"`
	TraceStatus          string    `json:"trace_status"`
	FailurePhase         string    `json:"failure_phase"`
	FailureReason        string    `json:"failure_reason"`
	RawEvidenceFile      string    `json:"raw_unverified_file"`
	RawEvidenceSHA256    string    `json:"raw_unverified_sha256"`
	PrimaryCommit        string    `json:"primary_commit"`
	FixturePrimaryCommit string    `json:"fixture_primary_commit"`
	ProductionFastCommit string    `json:"production_fast_commit"`
	DiagnosticHead       string    `json:"diagnostic_head"`
	CallBudget           int       `json:"call_budget"`
	ActualCalls          int       `json:"actual_calls"`
}

// TestFastDiagPublicE32Trace consumes the exact held Level0 ciphertext exported
// by Primary perfprobe. It executes one untraced cold call and one traced warm
// call only; the caller reserves both attempts before spawning this test.
func TestFastDiagPublicE32Trace(t *testing.T) {
	manifestPath := os.Getenv("FASTDIAG_PUBLIC_E32_MANIFEST")
	referenceVectorsPath := os.Getenv("FASTDIAG_PUBLIC_E32_FAST_VECTORS")
	outputPath := os.Getenv("FASTDIAG_OUTPUT")
	rawOutputPath := os.Getenv("FASTDIAG_RAW_OUTPUT")
	failureOutputPath := os.Getenv("FASTDIAG_FAILURE_OUTPUT")
	profileDir := os.Getenv("FASTDIAG_PROFILE_DIR")
	primaryHead := os.Getenv("FASTDIAG_PRIMARY_HEAD")
	scopes := os.Getenv("FASTDIAG_TRACE")
	require.NotEmpty(t, manifestPath)
	require.NotEmpty(t, referenceVectorsPath)
	require.NotEmpty(t, outputPath)
	require.NotEmpty(t, rawOutputPath)
	require.NotEmpty(t, failureOutputPath)
	require.NotEmpty(t, profileDir)
	require.NotEmpty(t, primaryHead)
	require.Equal(t, "stage,power,rescale", scopes, "E32 trace scope is frozen for this task")
	profileInfo, err := os.Stat(profileDir)
	require.NoError(t, err)
	require.True(t, profileInfo.IsDir())
	for _, path := range []string{outputPath, rawOutputPath, failureOutputPath, filepath.Join(profileDir, "warm-cpu.pprof"), filepath.Join(profileDir, "post-warm-heap.pprof")} {
		_, err := os.Lstat(path)
		require.True(t, os.IsNotExist(err), "refusing to overwrite diagnostic output %s (stat error: %v)", path, err)
	}

	manifestBytes, err := os.ReadFile(manifestPath)
	require.NoError(t, err)
	var manifest fastDiagPublicE32Manifest
	require.NoError(t, json.Unmarshal(manifestBytes, &manifest))
	require.Equal(t, "fast-public-e32-trace-fixture.v1", manifest.SchemaVersion)
	require.Equal(t, "logn13-e32-public-native", manifest.Profile)
	require.Equal(t, "fast", manifest.Backend)
	require.Equal(t, "2d6145d7e1db0ca7351eb47a03e1b352fc4ef9ac", manifest.BackendCommit)
	require.Equal(t, 32, manifest.EphemeralSecretWeight)
	require.Equal(t, 1e-6, manifest.NumericalGate)
	require.Equal(t, 0, manifest.State.Level)
	require.Equal(t, 1, manifest.State.Degree)
	require.Equal(t, 1, manifest.State.QPrefixRows)
	require.True(t, manifest.IsNTT)
	require.False(t, manifest.IsMontgomery)
	require.False(t, manifest.C1Nonzero)
	require.NotEmpty(t, manifest.CiphertextFingerprintSHA256)
	require.Len(t, manifest.PhysicalRowLengths, 2)

	fixtureDir := filepath.Dir(manifestPath)
	for _, name := range []string{manifest.ParametersFile, manifest.CiphertextFile} {
		require.NotEmpty(t, name)
		require.Equal(t, filepath.Base(name), name, "fixture files must be local basenames")
	}
	parameterBytes, err := os.ReadFile(filepath.Join(fixtureDir, manifest.ParametersFile))
	require.NoError(t, err)
	require.Equal(t, manifest.ParametersSHA256, fastDiagSHA256(parameterBytes))
	ciphertextBytes, err := os.ReadFile(filepath.Join(fixtureDir, manifest.CiphertextFile))
	require.NoError(t, err)
	require.Equal(t, manifest.CiphertextSHA256, fastDiagSHA256(ciphertextBytes))

	var params Parameters
	require.NoError(t, params.UnmarshalBinary(parameterBytes))
	require.Equal(t, 32, params.EphemeralSecretWeight)
	require.Equal(t, manifest.Parameters.LogN, params.BootstrappingParameters.LogN())
	require.Equal(t, manifest.Parameters.LogSlots, params.CoeffsToSlotsParameters.LogSlots)
	require.Equal(t, manifest.Parameters.RingN, params.BootstrappingParameters.N())
	require.Equal(t, ModUpThenEncode, params.CircuitOrder)
	require.Equal(t, 30, params.Mod1ParametersLiteral.Mod1Degree)
	require.Equal(t, 3, params.Mod1ParametersLiteral.DoubleAngle)
	require.Equal(t, 16, params.Mod1ParametersLiteral.K)
	require.Equal(t, 60, params.Mod1ParametersLiteral.LogScale)
	require.Equal(t, manifest.Parameters.QPrimes, fastDiagPrimeStrings(params.BootstrappingParameters.Q()))
	require.Equal(t, manifest.Parameters.PPrimes, fastDiagPrimeStrings(params.BootstrappingParameters.P()))
	require.Equal(t, strconv.FormatUint(params.BootstrappingParameters.Q()[0], 10), manifest.State.PrefixQ)
	require.Equal(t, 1<<12, manifest.Parameters.InputSlots)

	input := rlwe.NewCiphertext(params.ResidualParameters, 1, 0)
	require.NoError(t, input.UnmarshalBinary(ciphertextBytes))
	require.NotNil(t, input.MetaData)
	require.Equal(t, 0, input.Level())
	require.Equal(t, 1, input.Degree())
	require.Equal(t, manifest.State.Scale, input.Scale.Value.Text('e', 80))
	require.True(t, input.Scale.Equal(params.ResidualParameters.DefaultScale()))
	require.Equal(t, manifest.Parameters.LogSlots, input.LogSlots())
	require.Equal(t, manifest.IsNTT, input.IsNTT)
	require.Equal(t, manifest.IsMontgomery, input.IsMontgomery)
	rows, err := fastckks.QPrefixWidth(input.Level())
	require.NoError(t, err)
	require.Equal(t, 1, rows)
	for component := range input.Value {
		require.Len(t, input.Value[component].Coeffs, len(manifest.PhysicalRowLengths[component]))
		for row := range input.Value[component].Coeffs {
			require.Equal(t, manifest.PhysicalRowLengths[component][row], len(input.Value[component].Coeffs[row]))
		}
	}
	for _, residue := range input.Value[1].Coeffs[0] {
		require.Zero(t, residue, "Fast public trace input must remain in the accepted zero-secret mode")
	}
	require.Equal(t, manifest.CiphertextFingerprintSHA256, fastDiagPublicCiphertextFingerprint(t, input))
	decodedInput := fastDiagDecode(t, params.ResidualParameters, input)
	require.Equal(t, manifest.DecodedSHA256, fastDiagInputSHA256(decodedInput))

	vectorBytes, err := os.ReadFile(referenceVectorsPath)
	require.NoError(t, err)
	var vectors fastDiagPublicVectors
	require.NoError(t, json.Unmarshal(vectorBytes, &vectors))
	validateOnly := os.Getenv("FASTDIAG_VALIDATE_ONLY") == "1"
	if validateOnly {
		preflightVectors := vectors.SchemaVersion == "fast-standard-public-native-vectors.v1" && vectors.Mode == "public-native"
		repeatabilityVectors := vectors.SchemaVersion == "fast-standard-public-native-repeatability-vectors.v1" && vectors.Mode == "public-native-repeatability"
		require.True(t, preflightVectors || repeatabilityVectors, "validate-only mode requires a public-native held-input vector artifact")
	} else {
		require.Equal(t, "fast-standard-public-native-repeatability-vectors.v1", vectors.SchemaVersion)
		require.Equal(t, "public-native-repeatability", vectors.Mode)
	}
	require.Equal(t, "fast", vectors.Backend)
	require.Equal(t, manifest.BackendCommit, vectors.BackendCommit)
	require.Equal(t, manifest.PrimaryCommit, vectors.PrimaryCommit)
	require.Equal(t, manifest.PrimarySourceSHA256, vectors.SourceSHA256)
	require.Equal(t, manifest.ConfigSHA256, vectors.ConfigSHA256)
	require.Equal(t, manifest.QPSHA256, vectors.QPSHA256)
	require.Equal(t, manifest.InputSHA256, vectors.InputSHA256)
	require.Equal(t, manifest.WorkloadSHA256, vectors.WorkloadSHA256)
	inputOraclePairs, ok := vectors.Checkpoints["drop_level0"]
	require.True(t, ok)
	inputOracle := fastDiagDecodePairs(inputOraclePairs)
	require.Len(t, inputOracle, 1<<12)
	require.True(t, fastDiagFinite(inputOracle))
	require.Equal(t, manifest.DecodedSHA256, fastDiagInputSHA256(inputOracle))
	require.LessOrEqual(t, fastDiagVectorMetrics(inputOracle, decodedInput).max, manifest.NumericalGate)
	if validateOnly {
		// A zero-call preflight has only the held drop_level0 vector. Validate
		// its binary fixture compatibility without requiring Bootstrap outputs.
		return
	}
	fastReferenceNames := []string{"bootstrap_first_cold", "bootstrap_warm_01", "bootstrap_warm_02", "bootstrap_warm_03", "bootstrap_warm_04", "bootstrap_warm_05"}
	fastReferences := make([][]complex128, len(fastReferenceNames))
	for index, name := range fastReferenceNames {
		pairs, exists := vectors.BootstrapOutputs[name]
		require.True(t, exists, "missing uninstrumented Fast output %s", name)
		fastReferences[index] = fastDiagDecodePairs(pairs)
		require.Len(t, fastReferences[index], 1<<12)
		require.True(t, fastDiagFinite(fastReferences[index]))
	}
	require.Len(t, vectors.BootstrapOutputs, 6, "reference vectors must contain exactly six uninstrumented Fast outputs")

	keygen := rlwe.NewKeyGenerator(params.ResidualParameters)
	sk := keygen.GenSecretKeyNew()
	keys, _, err := params.GenEvaluationKeys(sk)
	require.NoError(t, err)
	require.True(t, keys.fastCompatible, "public GenEvaluationKeys must select Fast-compatible dispatch")
	eval, err := NewEvaluator(params, keys)
	require.NoError(t, err)
	require.NotNil(t, eval.fast, "ordinary NewEvaluator must dispatch to the intended Fast kernel")
	require.Nil(t, eval.Evaluator, "public Fast dispatch must not construct the Standard evaluator")

	inputBytesBefore, err := input.MarshalBinary()
	require.NoError(t, err)
	manifestSum := sha256.Sum256(manifestBytes)
	cpuProfilePath := filepath.Join(profileDir, "warm-cpu.pprof")
	heapProfilePath := filepath.Join(profileDir, "post-warm-heap.pprof")
	artifact := fastDiagPublicE32Artifact{
		SchemaVersion: "fastdiag.public-e32.trace.v1", Timestamp: time.Now().UTC(), TraceStatus: "TRACE_UNVERIFIED",
		TraceValidationError: "raw event capture is not certified; a separate certified result is written only after all gates pass",
		Profile:              manifest.Profile, Mode: "test-only-public-e32-fast", PrimaryCommit: primaryHead, FixturePrimaryCommit: manifest.PrimaryCommit,
		PrimarySourceSHA256: manifest.PrimarySourceSHA256, ProductionFastCommit: manifest.BackendCommit,
		DiagnosticHead: os.Getenv("FASTDIAG_DIAGNOSTIC_HEAD"), ConfigSHA256: manifest.ConfigSHA256,
		QPSHA256: manifest.QPSHA256, InputSHA256: manifest.InputSHA256, WorkloadSHA256: manifest.WorkloadSHA256,
		FixtureManifestSHA256: hex.EncodeToString(manifestSum[:]), CiphertextSHA256: manifest.CiphertextSHA256,
		ExpectedInputSHA256: manifest.DecodedSHA256, CallBudget: 2, NumericalGate: manifest.NumericalGate,
		GoVersion: runtime.Version(), OS: runtime.GOOS, Arch: runtime.GOARCH,
		NumCPU: runtime.NumCPU(), GOMAXPROCS: runtime.GOMAXPROCS(0),
		GOGC: fastDiagEnvironmentSetting("GOGC"), GOMEMLIMIT: fastDiagEnvironmentSetting("GOMEMLIMIT"), GODEBUG: fastDiagEnvironmentSetting("GODEBUG"),
		CPUProfileFile: filepath.Base(cpuProfilePath), HeapProfileFile: filepath.Base(heapProfilePath),
	}
	callCount := 0
	rawSHA256 := ""
	failureArtifactWritten := false
	writeRaw := func() error {
		artifact.Timestamp = time.Now().UTC()
		artifact.ActualCalls = callCount
		var writeErr error
		rawSHA256, writeErr = fastDiagWriteExclusiveJSON(rawOutputPath, artifact)
		return writeErr
	}
	writeFailure := func(phase, reason string) {
		if rawSHA256 == "" {
			if err := writeRaw(); err != nil {
				t.Errorf("cannot preserve unverified raw E32 trace %s: %v", rawOutputPath, err)
				return
			}
		}
		if err := fastDiagWritePublicE32Failure(failureOutputPath, rawOutputPath, rawSHA256, phase, reason, artifact); err != nil {
			t.Errorf("cannot preserve E32 trace failure reason at %s: %v", failureOutputPath, err)
			return
		}
		failureArtifactWritten = true
	}
	fastdiag.Reset()
	if err := fastdiag.ConfigureCSV(""); err != nil {
		t.Fatalf("configure untraced cold call: %v", err)
	}
	coldRun := fastDiagPublicE32Run{Index: 1, Phase: "first_cold_bootstrap"}
	coldStarted := time.Now()
	fastdiag.Reset()
	coldOutput, err := eval.Bootstrap(input.CopyNew())
	coldElapsed := time.Since(coldStarted).Nanoseconds()
	callCount++
	coldRun.ElapsedNS = coldElapsed
	coldOK := false
	if err == nil {
		coldOK = t.Run("cold_native_output_oracle", func(t *testing.T) {
			coldRun = fastDiagValidatePublicE32Output(t, params, manifest, coldOutput, inputOracle, fastReferences, 1, "first_cold_bootstrap", coldElapsed)
		})
	}
	coldEvents := fastdiag.Events()
	artifact.Runs = []fastDiagPublicE32Run{coldRun}
	artifact.Events = coldEvents
	if err != nil || !coldOK || len(coldEvents) != 0 {
		reason := "cold public Bootstrap output failed its native oracle or unexpectedly emitted trace events; warm call was not launched"
		if err != nil {
			reason = fmt.Sprintf("cold public Bootstrap failed: %v", err)
		}
		writeFailure("cold_call_oracle", reason)
		t.Fatalf("%s; raw trace: %s; failure record: %s", reason, rawOutputPath, failureOutputPath)
	}

	traceScopes := "stage,power,rescale"
	if err := fastdiag.ConfigureCSV(traceScopes); err != nil {
		writeFailure("trace_configuration", fmt.Sprintf("configure trace scopes: %v", err))
		t.Fatalf("configure trace scopes: %v; raw trace: %s", err, rawOutputPath)
	}
	fastdiag.Reset()
	cpuProfile, err := os.OpenFile(cpuProfilePath, os.O_WRONLY|os.O_CREATE|os.O_EXCL, 0o600)
	if err != nil {
		writeFailure("warm_profile_setup", fmt.Sprintf("create warm CPU profile: %v", err))
		t.Fatalf("create warm CPU profile: %v; raw trace: %s", err, rawOutputPath)
	}
	if err := pprof.StartCPUProfile(cpuProfile); err != nil {
		_ = cpuProfile.Close()
		writeFailure("warm_profile_setup", fmt.Sprintf("start warm CPU profile: %v", err))
		t.Fatalf("start warm CPU profile: %v; raw trace: %s", err, rawOutputPath)
	}
	warmStarted := time.Now()
	warmOutput, warmErr := eval.Bootstrap(input.CopyNew())
	warmElapsed := time.Since(warmStarted).Nanoseconds()
	pprof.StopCPUProfile()
	cpuCloseErr := cpuProfile.Close()
	callCount++
	warmRun := fastDiagPublicE32Run{Index: 2, Phase: "warm_bootstrap_01", ElapsedNS: warmElapsed}
	warmOracleOK := false
	if warmErr == nil {
		warmOracleOK = t.Run("warm_native_output_oracle", func(t *testing.T) {
			warmRun = fastDiagValidatePublicE32Output(t, params, manifest, warmOutput, inputOracle, fastReferences, 2, "warm_bootstrap_01", warmElapsed)
		})
	}
	events := fastdiag.Events()
	artifact.Runs = []fastDiagPublicE32Run{coldRun, warmRun}
	artifact.Events = events
	if err := writeRaw(); err != nil {
		t.Fatalf("cannot atomically preserve unverified E32 raw trace at %s: %v", rawOutputPath, err)
	}
	defer func() {
		if t.Failed() && !failureArtifactWritten {
			writeFailure("post_capture_validation", "a post-capture oracle, event, profile, fixture-integrity, or certification gate failed; see the test output")
		}
	}()
	if warmErr != nil {
		reason := fmt.Sprintf("warm public Bootstrap failed: %v", warmErr)
		writeFailure("warm_call", reason)
		t.Fatalf("%s; raw trace: %s; failure record: %s", reason, rawOutputPath, failureOutputPath)
	}
	if !warmOracleOK {
		reason := "warm decoded plaintext oracle or original-Fast reference gate failed"
		writeFailure("warm_output_oracle", reason)
		t.Fatalf("%s; raw trace: %s; failure record: %s", reason, rawOutputPath, failureOutputPath)
	}
	if cpuCloseErr != nil {
		reason := fmt.Sprintf("close warm CPU profile: %v", cpuCloseErr)
		writeFailure("warm_profile", reason)
		t.Fatalf("%s; raw trace: %s; failure record: %s", reason, rawOutputPath, failureOutputPath)
	}
	if eventErr := fastDiagPublicE32EventsValidationError(events); eventErr != nil {
		writeFailure("event_tree_validation", eventErr.Error())
		t.Fatalf("E32 event-tree validation failed: %v; raw TRACE_UNVERIFIED: %s; failure record: %s", eventErr, rawOutputPath, failureOutputPath)
	}

	heapProfile, err := os.OpenFile(heapProfilePath, os.O_WRONLY|os.O_CREATE|os.O_EXCL, 0o600)
	require.NoError(t, err)
	runtime.GC()
	heapWriteErr := pprof.WriteHeapProfile(heapProfile)
	heapCloseErr := heapProfile.Close()
	require.NoError(t, heapWriteErr)
	require.NoError(t, heapCloseErr)
	cpuProfileBytes, err := os.ReadFile(cpuProfilePath)
	require.NoError(t, err)
	heapProfileBytes, err := os.ReadFile(heapProfilePath)
	require.NoError(t, err)
	require.NotEmpty(t, cpuProfileBytes)
	require.NotEmpty(t, heapProfileBytes)
	require.Equal(t, 2, callCount, "the E32 diagnostic lane must execute exactly cold plus one traced warm call")
	inputBytesAfter, err := input.MarshalBinary()
	require.NoError(t, err)
	require.Equal(t, fastDiagSHA256(inputBytesBefore), fastDiagSHA256(inputBytesAfter), "Bootstrap must not mutate the held fixture")

	artifact.CPUProfileSHA256 = fastDiagSHA256(cpuProfileBytes)
	artifact.HeapProfileFile = filepath.Base(heapProfilePath)
	artifact.HeapProfileSHA256 = fastDiagSHA256(heapProfileBytes)
	artifact.TraceStatus = "TRACE_VALIDATED"
	artifact.TraceValidationError = ""
	artifact.RawEvidenceFile = filepath.Base(rawOutputPath)
	artifact.RawEvidenceSHA256 = rawSHA256
	artifact.Timestamp = time.Now().UTC()
	_, err = fastDiagWriteExclusiveJSON(outputPath, artifact)
	require.NoError(t, err)
}

func fastDiagEnvironmentSetting(name string) string {
	if value := os.Getenv(name); value != "" {
		return value
	}
	return "runtime-default"
}

func fastDiagWriteExclusiveJSON(path string, value any) (string, error) {
	data, err := json.MarshalIndent(value, "", "  ")
	if err != nil {
		return "", err
	}
	data = append(data, '\n')
	temp, err := os.CreateTemp(filepath.Dir(path), "."+filepath.Base(path)+".tmp-")
	if err != nil {
		return "", err
	}
	tempPath := temp.Name()
	defer os.Remove(tempPath)
	if err := temp.Chmod(0o600); err != nil {
		_ = temp.Close()
		return "", err
	}
	if _, err := temp.Write(data); err != nil {
		_ = temp.Close()
		return "", err
	}
	if err := temp.Sync(); err != nil {
		_ = temp.Close()
		return "", err
	}
	if err := temp.Close(); err != nil {
		return "", err
	}
	if err := os.Link(tempPath, path); err != nil {
		return "", err
	}
	// The destination already has the synced bytes and is owner-only. Cleanup
	// failure for the temporary hard link must not report failure after the
	// exclusive destination has been successfully published.
	_ = os.Remove(tempPath)
	return fastDiagSHA256(data), nil
}

func fastDiagWritePublicE32Failure(path, rawPath, rawSHA256, phase, reason string, artifact fastDiagPublicE32Artifact) error {
	if path == "" || rawPath == "" || rawSHA256 == "" || phase == "" || reason == "" {
		return errors.New("E32 diagnostic failure artifact requires paths, phase, and reason")
	}
	failure := fastDiagPublicE32FailureArtifact{
		SchemaVersion: "fastdiag.public-e32.failure.v1", Timestamp: time.Now().UTC(),
		TraceStatus: "TRACE_UNVERIFIED", FailurePhase: phase, FailureReason: reason,
		RawEvidenceFile: filepath.Base(rawPath), RawEvidenceSHA256: rawSHA256,
		PrimaryCommit: artifact.PrimaryCommit, FixturePrimaryCommit: artifact.FixturePrimaryCommit, ProductionFastCommit: artifact.ProductionFastCommit,
		DiagnosticHead: artifact.DiagnosticHead, CallBudget: artifact.CallBudget, ActualCalls: artifact.ActualCalls,
	}
	_, err := fastDiagWriteExclusiveJSON(path, failure)
	return err
}

func TestFastDiagWriteExclusiveJSONIsOwnerOnlyAndDoesNotOverwrite(t *testing.T) {
	path := filepath.Join(t.TempDir(), "trace.json")
	gotSHA, err := fastDiagWriteExclusiveJSON(path, map[string]string{"trace_status": "TRACE_UNVERIFIED"})
	require.NoError(t, err)
	data, err := os.ReadFile(path)
	require.NoError(t, err)
	require.Equal(t, fastDiagSHA256(data), gotSHA)
	info, err := os.Stat(path)
	require.NoError(t, err)
	require.Equal(t, os.FileMode(0o600), info.Mode().Perm())
	_, err = fastDiagWriteExclusiveJSON(path, map[string]string{"trace_status": "TRACE_VALIDATED"})
	require.Error(t, err, "exclusive evidence writer must not replace an existing raw artifact")
	dataAfter, err := os.ReadFile(path)
	require.NoError(t, err)
	require.Equal(t, data, dataAfter)
}

func TestFastDiagWritePublicE32FailureKeepsReasonSourceAndBudget(t *testing.T) {
	dir := t.TempDir()
	rawPath := filepath.Join(dir, "raw-trace-unverified.json")
	rawSHA, err := fastDiagWriteExclusiveJSON(rawPath, map[string]string{"trace_status": "TRACE_UNVERIFIED"})
	require.NoError(t, err)
	failurePath := filepath.Join(dir, "trace-failure.json")
	artifact := fastDiagPublicE32Artifact{
		PrimaryCommit: "primary-head", FixturePrimaryCommit: "fixture-primary",
		ProductionFastCommit: "production-fast", DiagnosticHead: "diagnostic-head",
		CallBudget: 2, ActualCalls: 2,
	}
	require.NoError(t, fastDiagWritePublicE32Failure(failurePath, rawPath, rawSHA, "event_tree_validation", "missing materialization span", artifact))
	data, err := os.ReadFile(failurePath)
	require.NoError(t, err)
	var failure fastDiagPublicE32FailureArtifact
	require.NoError(t, json.Unmarshal(data, &failure))
	require.Equal(t, "TRACE_UNVERIFIED", failure.TraceStatus)
	require.Equal(t, "event_tree_validation", failure.FailurePhase)
	require.Equal(t, "missing materialization span", failure.FailureReason)
	require.Equal(t, filepath.Base(rawPath), failure.RawEvidenceFile)
	require.Equal(t, rawSHA, failure.RawEvidenceSHA256)
	require.Equal(t, "primary-head", failure.PrimaryCommit)
	require.Equal(t, "fixture-primary", failure.FixturePrimaryCommit)
	require.Equal(t, 2, failure.CallBudget)
	require.Equal(t, 2, failure.ActualCalls)
}

func fastDiagValidatePublicE32Output(
	t testing.TB,
	params Parameters,
	manifest fastDiagPublicE32Manifest,
	output *rlwe.Ciphertext,
	oracle []complex128,
	references [][]complex128,
	index int,
	phase string,
	elapsedNS int64,
) fastDiagPublicE32Run {
	t.Helper()
	require.NotNil(t, output)
	require.NotNil(t, output.MetaData)
	require.Equal(t, params.ResidualParameters.MaxLevel(), output.Level())
	require.Equal(t, 1, output.Degree())
	require.True(t, output.Scale.Equal(params.ResidualParameters.DefaultScale()))
	require.True(t, output.IsNTT)
	require.False(t, output.IsMontgomery)
	rows, err := fastckks.QPrefixWidth(output.Level())
	require.NoError(t, err)
	decoded := fastDiagDecode(t, params.ResidualParameters, output)
	require.Len(t, decoded, 1<<12)
	require.True(t, fastDiagFinite(decoded))
	oracleMetrics := fastDiagVectorMetrics(oracle, decoded)
	require.LessOrEqual(t, oracleMetrics.max, manifest.NumericalGate, "%s plaintext oracle", phase)
	worstRMSE, worstMax := 0.0, 0.0
	for refIndex, reference := range references {
		metrics := fastDiagVectorMetrics(reference, decoded)
		if metrics.rmse > worstRMSE {
			worstRMSE = metrics.rmse
		}
		if metrics.max > worstMax {
			worstMax = metrics.max
		}
		require.LessOrEqual(t, metrics.max, manifest.NumericalGate, "%s vs uninstrumented Fast sample %d", phase, refIndex+1)
	}
	return fastDiagPublicE32Run{
		Index: index, Phase: phase, ElapsedNS: elapsedNS, DecodedSHA256: fastDiagInputSHA256(decoded),
		OracleRMSE: oracleMetrics.rmse, OracleMaxComplex: oracleMetrics.max,
		FastReferenceWorstRMSE: worstRMSE, FastReferenceWorstMax: worstMax,
		Level: output.Level(), Degree: output.Degree(), Scale: output.Scale.Value.Text('e', 80), QPrefixRows: rows,
	}
}

func fastDiagValidatePublicE32Events(t testing.TB, events []fastdiag.Event) {
	t.Helper()
	require.NoError(t, fastDiagPublicE32EventsValidationError(events))
}

func fastDiagPublicE32EventsValidationError(events []fastdiag.Event) error {
	if len(events) == 0 {
		return errors.New("traced warm call contains no fastdiag events")
	}
	rootSequence := uint64(0)
	var stages []fastdiag.Event
	generatedPowerParents := map[uint64]bool{}
	children := map[uint64][]fastdiag.Event{}
	sequences := map[uint64]bool{}
	for _, event := range events {
		if event.Sequence == 0 || event.ElapsedNS < 0 {
			return fmt.Errorf("fastdiag event %q/%q has zero sequence or negative elapsed time", event.Scope, event.Name)
		}
		if sequences[event.Sequence] {
			return fmt.Errorf("duplicate fastdiag event sequence %d", event.Sequence)
		}
		sequences[event.Sequence] = true
		children[event.ParentSequence] = append(children[event.ParentSequence], event)
		if event.Scope == fastdiag.Stage {
			if event.Name == "bootstrap" {
				if rootSequence != 0 || event.ParentSequence != 0 {
					return errors.New("exactly one root Bootstrap event is expected")
				}
				rootSequence = event.Sequence
			} else {
				stages = append(stages, event)
			}
		}
		if event.Scope == fastdiag.Power && event.Name == "generated_powers" {
			generatedPowerParents[event.Sequence] = true
		}
	}
	if rootSequence == 0 {
		return errors.New("missing Bootstrap root event")
	}
	for _, event := range events {
		if event.ParentSequence != 0 {
			if !sequences[event.ParentSequence] {
				return fmt.Errorf("event %d references missing parent %d", event.Sequence, event.ParentSequence)
			}
		}
	}
	stageNames := make([]string, len(stages))
	for index, stage := range stages {
		if stage.ParentSequence != rootSequence {
			return fmt.Errorf("stage %s must be a direct child of Bootstrap", stage.Name)
		}
		stageNames[index] = stage.Name
	}
	wantStages := []string{"pack_n1_to_n2", "scale_down", "mod_up_trace", "coeffs_to_slots", "evalmod_real"}
	for _, optional := range []string{"evalmod_imag"} {
		found := false
		for _, name := range stageNames {
			found = found || name == optional
		}
		if found {
			wantStages = append(wantStages, optional)
		}
	}
	wantStages = append(wantStages, "slots_to_coeffs", "unpack_n2_to_n1", "public_finalization")
	if len(wantStages) != len(stageNames) {
		return fmt.Errorf("stage sequence mismatch: got %v, want %v", stageNames, wantStages)
	}
	for index := range wantStages {
		if stageNames[index] != wantStages[index] {
			return fmt.Errorf("stage sequence mismatch: got %v, want %v", stageNames, wantStages)
		}
	}
	powerEvents, rescaleParents := 0, 0
	generatedPowers := map[int]bool{}
	for _, event := range events {
		switch event.Scope {
		case fastdiag.Power:
			if event.Name == "power" {
				powerEvents++
				if !generatedPowerParents[event.ParentSequence] {
					return errors.New("power event is not nested under generated_powers")
				}
				if event.Power == nil || event.SplitA == nil || event.SplitB == nil {
					return errors.New("generated power event is missing its power or recurrence split")
				}
				generatedPowers[*event.Power] = true
			}
		case fastdiag.Rescale:
			if event.Name == "rescale" {
				rescaleParents++
				passChildren := children[event.Sequence]
				if len(passChildren) != 2 {
					return fmt.Errorf("Rescale event %d must include exactly preflight and materialization children", event.Sequence)
				}
				passes := make(map[string]fastdiag.Event, 2)
				for _, pass := range passChildren {
					if pass.Scope != fastdiag.Rescale {
						return fmt.Errorf("Rescale event %d has a non-Rescale direct child", event.Sequence)
					}
					if _, exists := passes[pass.Name]; exists {
						return fmt.Errorf("Rescale event %d has duplicate child %q", event.Sequence, pass.Name)
					}
					passes[pass.Name] = pass
				}
				preflight, ok := passes["preflight"]
				if !ok {
					return fmt.Errorf("Rescale event %d is missing preflight child", event.Sequence)
				}
				materialization, ok := passes["materialization"]
				if !ok {
					return fmt.Errorf("Rescale event %d is missing materialization child", event.Sequence)
				}
				if err := fastDiagPublicE32RescalePassError(preflight, children, false); err != nil {
					return err
				}
				if err := fastDiagPublicE32RescalePassError(materialization, children, true); err != nil {
					return err
				}
			}
		}
	}
	if powerEvents == 0 {
		return errors.New("trace has no generated power events")
	}
	if rescaleParents == 0 {
		return errors.New("trace has no Rescale parent events")
	}
	for _, power := range []int{2, 3, 4, 6, 8, 16} {
		if !generatedPowers[power] {
			return fmt.Errorf("missing generated Chebyshev power T%d", power)
		}
	}
	return nil
}

func fastDiagPublicE32RescalePassError(pass fastdiag.Event, children map[uint64][]fastdiag.Event, materializing bool) error {
	if pass.Count == nil {
		return fmt.Errorf("Rescale %s pass is missing its component count", pass.Name)
	}
	expectedComponents := *pass.Count
	passChildren := children[pass.Sequence]
	prefixByComponent := map[string]int{}
	restoreByComponent := map[string]int{}
	loopsByComponent := map[string]fastdiag.Event{}
	for _, child := range passChildren {
		if child.Scope != fastdiag.Rescale {
			return fmt.Errorf("Rescale %s has non-Rescale child %q", pass.Name, child.Name)
		}
		switch child.Name {
		case "prefix_to_coefficient":
			prefixByComponent[child.Component]++
		case "coefficient_loop":
			if _, exists := loopsByComponent[child.Component]; exists {
				return fmt.Errorf("duplicate coefficient loop for %s/%s", pass.Name, child.Component)
			}
			loopsByComponent[child.Component] = child
		case "ntt_montgomery_restore":
			if !materializing {
				return errors.New("preflight must not contain an NTT restore event")
			}
			restoreByComponent[child.Component]++
		default:
			return fmt.Errorf("unexpected direct child %q under Rescale %s", child.Name, pass.Name)
		}
	}
	if materializing {
		if len(prefixByComponent) != 0 || len(loopsByComponent) != 0 || len(restoreByComponent) != expectedComponents {
			return fmt.Errorf("Rescale materialization must contain exactly %d per-component restore events", expectedComponents)
		}
		for component, count := range restoreByComponent {
			if component == "" || count != 1 {
				return fmt.Errorf("Rescale materialization has invalid restore count for component %q", component)
			}
		}
	} else {
		if len(prefixByComponent) != expectedComponents || len(loopsByComponent) != expectedComponents || len(restoreByComponent) != 0 {
			return fmt.Errorf("Rescale preflight must contain prefix and coefficient-loop events for %d components only", expectedComponents)
		}
	}
	for component, count := range prefixByComponent {
		if count != 1 {
			return fmt.Errorf("duplicate prefix conversion for %s/%s", pass.Name, component)
		}
		loop, ok := loopsByComponent[component]
		if !ok {
			return fmt.Errorf("Rescale %s component %s is missing its coefficient loop", pass.Name, component)
		}
		loopNames := map[string]int{}
		for _, child := range children[loop.Sequence] {
			loopNames[child.Name]++
		}
		if len(loopNames) != 2 || loopNames["reconstruct_center_round_capacity"] != 1 || loopNames["residue_materialization"] != 1 {
			return fmt.Errorf("Rescale %s component %s has an invalid coefficient-loop child set", pass.Name, component)
		}
	}
	return nil
}

func TestFastDiagPublicE32EventValidatorAcceptsSourceShapedTree(t *testing.T) {
	require.NoError(t, fastDiagPublicE32EventsValidationError(fastDiagPublicE32SourceShapedTree()))
}

func TestFastDiagPublicE32EventValidatorReportsMissingMaterialization(t *testing.T) {
	events := fastDiagPublicE32SourceShapedTree()
	var materializationSequence uint64
	for _, event := range events {
		if event.Name == "materialization" {
			materializationSequence = event.Sequence
		}
	}
	require.NotZero(t, materializationSequence)
	filtered := make([]fastdiag.Event, 0, len(events))
	for _, event := range events {
		if event.Sequence != materializationSequence && event.ParentSequence != materializationSequence {
			filtered = append(filtered, event)
		}
	}
	err := fastDiagPublicE32EventsValidationError(filtered)
	require.ErrorContains(t, err, "must include exactly preflight and materialization children")
}

func fastDiagPublicE32SourceShapedTree() []fastdiag.Event {
	stageNames := []string{"bootstrap", "pack_n1_to_n2", "scale_down", "mod_up_trace", "coeffs_to_slots", "evalmod_real", "evalmod_imag", "slots_to_coeffs", "unpack_n2_to_n1", "public_finalization"}
	events := make([]fastdiag.Event, 0, 32)
	sequence := uint64(1)
	for _, name := range stageNames {
		parent := uint64(0)
		if name != "bootstrap" {
			parent = 1
		}
		events = append(events, fastdiag.Event{Scope: fastdiag.Stage, Name: name, Sequence: sequence, ParentSequence: parent})
		sequence++
	}
	events = append(events, fastdiag.Event{Scope: fastdiag.Power, Name: "generated_powers", Sequence: sequence})
	generatedPowersSequence := sequence
	sequence++
	for _, power := range []int{2, 3, 4, 6, 8, 16} {
		p, splitA, splitB := power, power/2, power/2
		events = append(events, fastdiag.Event{
			Scope: fastdiag.Power, Name: "power", Sequence: sequence, ParentSequence: generatedPowersSequence,
			Fields: fastdiag.Fields{Power: &p, SplitA: &splitA, SplitB: &splitB},
		})
		sequence++
	}
	rescaleSequence := sequence
	events = append(events, fastdiag.Event{Scope: fastdiag.Rescale, Name: "rescale", Sequence: rescaleSequence, ParentSequence: 6})
	sequence++
	preflightSequence, materializationSequence := sequence, sequence+1
	preflightCount, materializationCount := 2, 2
	events = append(events,
		fastdiag.Event{Scope: fastdiag.Rescale, Name: "preflight", Sequence: preflightSequence, ParentSequence: rescaleSequence, Fields: fastdiag.Fields{Count: &preflightCount}},
		fastdiag.Event{Scope: fastdiag.Rescale, Name: "materialization", Sequence: materializationSequence, ParentSequence: rescaleSequence, Fields: fastdiag.Fields{Count: &materializationCount}},
	)
	sequence += 2
	for _, component := range []string{"c0", "c1"} {
		componentName := component
		prefixSequence, loopSequence := sequence, sequence+1
		events = append(events,
			fastdiag.Event{Scope: fastdiag.Rescale, Name: "prefix_to_coefficient", Sequence: prefixSequence, ParentSequence: preflightSequence, Fields: fastdiag.Fields{Component: componentName}},
			fastdiag.Event{Scope: fastdiag.Rescale, Name: "coefficient_loop", Sequence: loopSequence, ParentSequence: preflightSequence, Fields: fastdiag.Fields{Component: componentName}},
		)
		sequence += 2
		events = append(events,
			fastdiag.Event{Scope: fastdiag.Rescale, Name: "reconstruct_center_round_capacity", Sequence: sequence, ParentSequence: loopSequence},
			fastdiag.Event{Scope: fastdiag.Rescale, Name: "residue_materialization", Sequence: sequence + 1, ParentSequence: loopSequence},
			fastdiag.Event{Scope: fastdiag.Rescale, Name: "ntt_montgomery_restore", Sequence: sequence + 2, ParentSequence: materializationSequence, Fields: fastdiag.Fields{Component: componentName}},
		)
		sequence += 3
	}
	return events
}

func fastDiagDecodePairs(pairs []fastDiagPublicVector) []complex128 {
	values := make([]complex128, len(pairs))
	for index, pair := range pairs {
		values[index] = complex(pair.Real, pair.Imag)
	}
	return values
}

func fastDiagVectorMetrics(want, got []complex128) fastDiagPublicMetrics {
	if len(want) == 0 || len(want) != len(got) {
		return fastDiagPublicMetrics{rmse: math.Inf(1), max: math.Inf(1)}
	}
	var sum float64
	var maxDifference float64
	for index := range want {
		difference := cmplxAbs(want[index] - got[index])
		sum += difference * difference
		maxDifference = math.Max(maxDifference, difference)
	}
	return fastDiagPublicMetrics{rmse: math.Sqrt(sum / float64(len(want))), max: maxDifference}
}

func fastDiagFinite(values []complex128) bool {
	for _, value := range values {
		if math.IsNaN(real(value)) || math.IsNaN(imag(value)) || math.IsInf(real(value), 0) || math.IsInf(imag(value), 0) {
			return false
		}
	}
	return len(values) != 0
}

func fastDiagSHA256(data []byte) string {
	sum := sha256.Sum256(data)
	return hex.EncodeToString(sum[:])
}

func fastDiagPublicCiphertextFingerprint(t testing.TB, ct *rlwe.Ciphertext) string {
	t.Helper()
	require.NotNil(t, ct)
	require.NotNil(t, ct.MetaData)
	hash := sha256.New()
	var encoded [8]byte
	writeUint64 := func(value uint64) {
		binary.LittleEndian.PutUint64(encoded[:], value)
		_, _ = hash.Write(encoded[:])
	}
	writeString := func(value string) {
		writeUint64(uint64(len(value)))
		_, _ = hash.Write([]byte(value))
	}
	writeUint64(uint64(ct.Level()))
	writeUint64(uint64(ct.Degree()))
	writeString(ct.Scale.Value.Text('e', 80))
	writeUint64(uint64(ct.LogDimensions.Rows))
	writeUint64(uint64(ct.LogDimensions.Cols))
	if ct.IsNTT {
		writeUint64(1)
	} else {
		writeUint64(0)
	}
	if ct.IsMontgomery {
		writeUint64(1)
	} else {
		writeUint64(0)
	}
	writeUint64(uint64(len(ct.Value)))
	for component, polynomial := range ct.Value {
		writeUint64(uint64(component))
		writeUint64(uint64(len(polynomial.Coeffs)))
		for row, coefficients := range polynomial.Coeffs {
			writeUint64(uint64(row))
			writeUint64(uint64(len(coefficients)))
			for _, coefficient := range coefficients {
				writeUint64(coefficient)
			}
		}
	}
	return hex.EncodeToString(hash.Sum(nil))
}

func fastDiagPrimeStrings(primes []uint64) []string {
	values := make([]string, len(primes))
	for index, prime := range primes {
		values[index] = strconv.FormatUint(prime, 10)
	}
	return values
}

func cmplxAbs(value complex128) float64 { return math.Hypot(real(value), imag(value)) }
