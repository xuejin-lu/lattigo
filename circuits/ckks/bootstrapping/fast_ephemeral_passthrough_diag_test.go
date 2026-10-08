//go:build fast_ephemeral_diag

package bootstrapping

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"math"
	"math/bits"
	"os"
	"runtime"
	"strconv"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const (
	fastEPassStandardCommit = "5dbffbdea05394de2ca3a432ed5318aa832e3f40"
	fastEPassConfigSHA256   = "919a2d9409b8ddeb720458d3112855f87d7a69e0cc39aad825b0ddf794769c98"
	fastEPassInputSHA256    = "d9151964e398ae9fb77248394b4b28f84c9e5737c621cf5e5b0570e343dcc285"
	fastEPassEffectiveSHA   = "f151442a4e08e1ebf8b7bb515fdd075a3748298bee1f3e2b990a6ad99aeca693"
	fastEPassQPrimesSHA     = "30066778ef1caea959b3357b0a70ff792fe8582ba9545f6b53c8df2a1a6cb568"
	fastEPassPPrimesSHA     = "70451211d27cd1e6bf092aa3c0752a63147c34842211dcd70d64d50755d318f7"
)

type fastEPassComplexPair struct {
	Real float64 `json:"real"`
	Imag float64 `json:"imag"`
}

type fastEPassMetrics struct {
	ComplexRMSE     float64 `json:"complex_rmse"`
	MaxComplexError float64 `json:"max_complex_error"`
	MaxRealError    float64 `json:"max_real_error"`
	MaxImagError    float64 `json:"max_imag_error"`
	SNRdB           float64 `json:"snr_db"`
}

type fastEPassMetadata struct {
	Level          int     `json:"level"`
	ScaleLog2      float64 `json:"scale_log2"`
	Degree         int     `json:"degree"`
	N              int     `json:"n"`
	ComponentCount int     `json:"component_count"`
	Slots          int     `json:"slots"`
	IsNTT          bool    `json:"is_ntt"`
	IsMontgomery   bool    `json:"is_montgomery"`
}

type fastEPassReference struct {
	SchemaVersion             string                 `json:"schema_version"`
	StandardCommit            string                 `json:"standard_commit"`
	StandardDirty             bool                   `json:"standard_dirty"`
	ConfigSHA256              string                 `json:"config_sha256"`
	InputSHA256               string                 `json:"input_sha256"`
	EffectiveParametersSHA256 string                 `json:"effective_parameters_sha256"`
	QPrimesSHA256             string                 `json:"q_primes_sha256"`
	PPrimesSHA256             string                 `json:"p_primes_sha256"`
	OutputMetadata            fastEPassMetadata      `json:"output_metadata"`
	OutputAgainstOriginal     fastEPassMetrics       `json:"output_against_original"`
	Values                    []fastEPassComplexPair `json:"values"`
}

type fastEPassEffectiveParameters struct {
	LogN             int      `json:"log_n"`
	LogSlots         int      `json:"log_slots"`
	InputSlots       int      `json:"input_slots"`
	RingN            int      `json:"ring_n"`
	ResidualRingN    int      `json:"residual_ring_n"`
	Q0Target         int      `json:"q0_config_target_bits"`
	Q0Bits           int      `json:"q0_bits"`
	QChainBits       []int    `json:"q_chain_bits"`
	PBits            []int    `json:"p_bits"`
	QPrimes          []string `json:"q_primes"`
	PPrimes          []string `json:"p_primes"`
	DefaultScale     string   `json:"default_scale"`
	Mod1Scale        int      `json:"mod1_log_scale"`
	Mod1Degree       int      `json:"mod1_degree"`
	DoubleAngle      int      `json:"double_angle"`
	K                int      `json:"k"`
	LogMessageRatio  int      `json:"log_message_ratio"`
	CircuitOrder     int      `json:"circuit_order"`
	QPrefixRowsAtMax int      `json:"q_prefix_rows_at_max_level"`
}

type fastEPassRun struct {
	Weight                  int               `json:"ephemeral_secret_weight"`
	DenseToSparseKeyPresent bool              `json:"dense_to_sparse_key_present"`
	SparseToDenseKeyPresent bool              `json:"sparse_to_dense_key_present"`
	PublicBootstrapCalls    int               `json:"public_bootstrap_calls"`
	OutputMetadata          fastEPassMetadata `json:"output_metadata"`
	OutputAgainstOriginal   fastEPassMetrics  `json:"output_against_original"`
}

type fastEPassEvidence struct {
	SchemaVersion             string           `json:"schema_version"`
	Task                      string           `json:"task"`
	SecondaryCommit           string           `json:"secondary_commit"`
	GoVersion                 string           `json:"go_version"`
	OS                        string           `json:"os"`
	Arch                      string           `json:"arch"`
	ConfigSHA256              string           `json:"config_sha256"`
	InputSHA256               string           `json:"input_sha256"`
	EffectiveParametersSHA256 string           `json:"effective_parameters_sha256"`
	QPrimesSHA256             string           `json:"q_primes_sha256"`
	PPrimesSHA256             string           `json:"p_primes_sha256"`
	InputKind                 string           `json:"input_kind"`
	Runs                      []fastEPassRun   `json:"runs"`
	E32VsE0                   fastEPassMetrics `json:"e32_vs_e0"`
	FastE32VsStandardE32      fastEPassMetrics `json:"fast_e32_vs_standard_e32"`
	StandardVsOriginal        fastEPassMetrics `json:"standard_e32_vs_original"`
}

// TestFastEphemeralPassthroughP93 runs exactly one public Fast Bootstrap at
// each E value. The Standard E=32 output is produced by a separate clean,
// pinned checkout and passed in as an aggregate-only comparison reference.
func TestFastEphemeralPassthroughP93(t *testing.T) {
	backendCommit := os.Getenv("FAST_E_PASSTHROUGH_BACKEND_COMMIT")
	require.Len(t, backendCommit, 40, "set the committed Secondary SHA used for this run")
	configPath := os.Getenv("FAST_E_PASSTHROUGH_CONFIG")
	require.NotEmpty(t, configPath)
	configBytes, err := os.ReadFile(configPath)
	require.NoError(t, err)
	configSHA := fastEPassHash(configBytes)
	require.Equal(t, fastEPassConfigSHA256, configSHA)

	standardReferencePath := os.Getenv("FAST_E_PASSTHROUGH_STANDARD_REFERENCE")
	require.NotEmpty(t, standardReferencePath)
	standardBytes, err := os.ReadFile(standardReferencePath)
	require.NoError(t, err)
	var standard fastEPassReference
	require.NoError(t, json.Unmarshal(standardBytes, &standard))
	require.Equal(t, "fast-zero-secret-e-passthrough-standard-reference.v1", standard.SchemaVersion)
	require.Equal(t, fastEPassStandardCommit, standard.StandardCommit)
	require.False(t, standard.StandardDirty)
	require.Equal(t, configSHA, standard.ConfigSHA256)

	params, residual, values, input := fastDiagP93Fixture(t)
	inputSHA := fastEPassInputHash(values)
	require.Equal(t, fastEPassInputSHA256, inputSHA)
	effective := fastEPassEffective(params, residual)
	effectiveJSON, err := json.Marshal(effective)
	require.NoError(t, err)
	effectiveSHA := fastEPassHash(effectiveJSON)
	require.Equal(t, fastEPassEffectiveSHA, effectiveSHA)
	qPrimesJSON, err := json.Marshal(effective.QPrimes)
	require.NoError(t, err)
	pPrimesJSON, err := json.Marshal(effective.PPrimes)
	require.NoError(t, err)
	qPrimesSHA, pPrimesSHA := fastEPassHash(qPrimesJSON), fastEPassHash(pPrimesJSON)
	require.Equal(t, fastEPassQPrimesSHA, qPrimesSHA)
	require.Equal(t, fastEPassPPrimesSHA, pPrimesSHA)
	require.Equal(t, effectiveSHA, standard.EffectiveParametersSHA256)
	require.Equal(t, qPrimesSHA, standard.QPrimesSHA256)
	require.Equal(t, pPrimesSHA, standard.PPrimesSHA256)
	require.Equal(t, inputSHA, standard.InputSHA256)
	require.Len(t, standard.Values, len(values))

	outputs := make(map[int][]complex128, 2)
	runs := make([]fastEPassRun, 0, 2)
	for _, weight := range []int{0, 32} {
		diagnosticParams := params
		diagnosticParams.EphemeralSecretWeight = weight
		sk := rlwe.NewKeyGenerator(residual).GenSecretKeyNew()
		keys, _, keyErr := diagnosticParams.GenEvaluationKeys(sk)
		require.NoError(t, keyErr)
		require.True(t, keys.fastCompatible, "the public GenEvaluationKeys API must select Fast-compatible key material")
		densePresent, sparsePresent := keys.EvkDenseToSparse != nil, keys.EvkSparseToDense != nil
		require.Equal(t, weight == 32, densePresent)
		require.Equal(t, weight == 32, sparsePresent)

		eval, evalErr := NewEvaluator(diagnosticParams, keys)
		require.NoError(t, evalErr)
		require.NotNil(t, eval.fast)
		require.Nil(t, eval.Evaluator, "the public Fast path must not construct a Standard evaluator")
		require.Equal(t, weight, eval.Parameters.EphemeralSecretWeight)
		require.Equal(t, weight, eval.fast.Parameters.EphemeralSecretWeight)

		calls := 0
		calls++
		out, bootstrapErr := eval.Bootstrap(input.CopyNew())
		require.NoError(t, bootstrapErr)
		require.NotNil(t, out)
		require.Equal(t, 1, calls, "exactly one public Fast Bootstrap is permitted per E value")
		decoded := fastEPassDecode(t, residual, out)
		require.Len(t, decoded, len(values))
		require.True(t, fastEPassFiniteValues(decoded))
		outputs[weight] = decoded
		runs = append(runs, fastEPassRun{
			Weight: weight, DenseToSparseKeyPresent: densePresent, SparseToDenseKeyPresent: sparsePresent,
			PublicBootstrapCalls: calls, OutputMetadata: fastEPassOutputMetadata(out, len(decoded)),
			OutputAgainstOriginal: fastEPassCompare(values, decoded),
		})
	}
	require.Len(t, runs, 2)
	require.Equal(t, 0, runs[0].Weight)
	require.Equal(t, 32, runs[1].Weight)
	for _, run := range runs {
		require.Equal(t, 1, run.PublicBootstrapCalls)
	}

	standardValues := make([]complex128, len(standard.Values))
	for i, value := range standard.Values {
		standardValues[i] = complex(value.Real, value.Imag)
	}
	require.True(t, fastEPassFiniteValues(standardValues))
	standardMetrics := fastEPassCompare(values, standardValues)
	fastStandardMetrics := fastEPassCompare(standardValues, outputs[32])
	require.InDelta(t, standardMetrics.ComplexRMSE, standard.OutputAgainstOriginal.ComplexRMSE, 1e-15)

	evidence := fastEPassEvidence{
		SchemaVersion: "fast-zero-secret-e-passthrough.v1", Task: "FAST-ZERO-SECRET-E-PASSTHROUGH-001",
		SecondaryCommit: backendCommit, GoVersion: runtime.Version(), OS: runtime.GOOS, Arch: runtime.GOARCH,
		ConfigSHA256: configSHA, InputSHA256: inputSHA, EffectiveParametersSHA256: effectiveSHA,
		QPrimesSHA256: qPrimesSHA, PPrimesSHA256: pPrimesSHA, InputKind: "diagnostic_plaintext_like_c1_zero",
		Runs: runs, E32VsE0: fastEPassCompare(outputs[0], outputs[32]),
		FastE32VsStandardE32: fastStandardMetrics, StandardVsOriginal: standardMetrics,
	}
	encoded, err := json.Marshal(evidence)
	require.NoError(t, err)
	t.Log(string(encoded))
}

func fastEPassEffective(params Parameters, residual ckks.Parameters) fastEPassEffectiveParameters {
	bootstrap := params.BootstrappingParameters
	q, p := bootstrap.Q(), bootstrap.P()
	qStrings, pStrings := fastEPassPrimeStrings(q), fastEPassPrimeStrings(p)
	qBits, pBits := fastEPassPrimeBits(q), fastEPassPrimeBits(p)
	logSlots := params.CoeffsToSlotsParameters.LogSlots
	defaultScale := residual.DefaultScale()
	return fastEPassEffectiveParameters{
		LogN: bootstrap.LogN(), LogSlots: logSlots, InputSlots: 1 << logSlots,
		RingN: bootstrap.N(), ResidualRingN: residual.N(), Q0Target: 55, Q0Bits: bits.Len64(q[0]),
		QChainBits: qBits, PBits: pBits, QPrimes: qStrings, PPrimes: pStrings,
		DefaultScale: defaultScale.Value.Text('e', 80), Mod1Scale: 60, Mod1Degree: 30,
		DoubleAngle: 3, K: 16, LogMessageRatio: 10, CircuitOrder: int(params.CircuitOrder),
		QPrefixRowsAtMax: min(bootstrap.MaxLevel()+1, 4),
	}
}

func fastEPassPrimeStrings(primes []uint64) []string {
	out := make([]string, len(primes))
	for i, prime := range primes {
		out[i] = strconv.FormatUint(prime, 10)
	}
	return out
}

func fastEPassPrimeBits(primes []uint64) []int {
	out := make([]int, len(primes))
	for i, prime := range primes {
		out[i] = bits.Len64(prime)
	}
	return out
}

func fastEPassInputHash(values []complex128) string {
	hash := sha256.New()
	var encoded [16]byte
	for _, value := range values {
		binary.LittleEndian.PutUint64(encoded[:8], math.Float64bits(real(value)))
		binary.LittleEndian.PutUint64(encoded[8:], math.Float64bits(imag(value)))
		_, _ = hash.Write(encoded[:])
	}
	return hex.EncodeToString(hash.Sum(nil))
}

func fastEPassHash(data []byte) string {
	sum := sha256.Sum256(data)
	return hex.EncodeToString(sum[:])
}

func fastEPassDecode(t *testing.T, params ckks.Parameters, ct *rlwe.Ciphertext) []complex128 {
	t.Helper()
	plaintext := ckks.NewPlaintext(params, ct.Level())
	*plaintext.MetaData = *ct.MetaData
	plaintext.Value.Copy(ct.Value[0])
	plaintext.IsNTT, plaintext.IsMontgomery = ct.IsNTT, ct.IsMontgomery
	values := make([]complex128, 1<<ct.LogSlots())
	require.NoError(t, ckks.NewEncoder(params).Decode(plaintext, values))
	return values
}

func fastEPassOutputMetadata(ct *rlwe.Ciphertext, slots int) fastEPassMetadata {
	return fastEPassMetadata{
		Level: ct.Level(), ScaleLog2: ct.Scale.Log2(), Degree: ct.Degree(), N: ct.N(),
		ComponentCount: len(ct.Value), Slots: slots, IsNTT: ct.IsNTT, IsMontgomery: ct.IsMontgomery,
	}
}

func fastEPassCompare(reference, actual []complex128) fastEPassMetrics {
	var squaredError, signalEnergy float64
	var maxComplex, maxReal, maxImag float64
	for i, want := range reference {
		deltaReal := real(actual[i] - want)
		deltaImag := imag(actual[i] - want)
		realError, imagError := math.Abs(deltaReal), math.Abs(deltaImag)
		complexError := math.Hypot(deltaReal, deltaImag)
		squaredError += complexError * complexError
		signalEnergy += real(want)*real(want) + imag(want)*imag(want)
		maxComplex = math.Max(maxComplex, complexError)
		maxReal = math.Max(maxReal, realError)
		maxImag = math.Max(maxImag, imagError)
	}
	count := float64(len(reference))
	rmse := math.Sqrt(squaredError / count)
	snr := math.Inf(1)
	if squaredError > 0 {
		snr = 10 * math.Log10(signalEnergy/squaredError)
	}
	return fastEPassMetrics{ComplexRMSE: rmse, MaxComplexError: maxComplex, MaxRealError: maxReal, MaxImagError: maxImag, SNRdB: snr}
}

func fastEPassFiniteValues(values []complex128) bool {
	for _, value := range values {
		if math.IsNaN(real(value)) || math.IsInf(real(value), 0) || math.IsNaN(imag(value)) || math.IsInf(imag(value), 0) {
			return false
		}
	}
	return true
}
