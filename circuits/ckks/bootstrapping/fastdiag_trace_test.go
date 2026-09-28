//go:build fastdiag

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
	"strings"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/internal/fastdiag"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const fastDiagAcceptedInputSHA256 = "d9151964e398ae9fb77248394b4b28f84c9e5737c621cf5e5b0570e343dcc285"

type fastDiagRun struct {
	Index          int              `json:"index"`
	ElapsedNS      int64            `json:"elapsed_ns"`
	NumericalMatch bool             `json:"numerical_match"`
	Events         []fastdiag.Event `json:"events"`
}

type fastDiagTraceArtifact struct {
	SchemaVersion string        `json:"schema_version"`
	Timestamp     time.Time     `json:"timestamp"`
	Profile       string        `json:"profile"`
	Trace         []string      `json:"trace"`
	Warmup        int           `json:"warmup"`
	Repetitions   int           `json:"repetitions"`
	InputSHA256   string        `json:"input_sha256"`
	InputLength   int           `json:"input_length"`
	LogN          int           `json:"log_n"`
	LogSlots      int           `json:"log_slots"`
	QChainBits    []int         `json:"q_chain_bits"`
	PBits         []int         `json:"p_bits"`
	Polynomial    int           `json:"polynomial_degree"`
	DoubleAngle   int           `json:"double_angle"`
	GoVersion     string        `json:"go_version"`
	OS            string        `json:"os"`
	Arch          string        `json:"arch"`
	NumCPU        int           `json:"num_cpu"`
	GOMAXPROCS    int           `json:"gomaxprocs"`
	Runs          []fastDiagRun `json:"runs"`
}

func TestFastDiagP93Q55Trace(t *testing.T) {
	scopes := os.Getenv("FASTDIAG_TRACE")
	if scopes == "" {
		scopes = "stage,power,rescale"
	}
	require.NoError(t, fastdiag.ConfigureCSV(scopes))
	warmup := fastDiagEnvInt(t, "FASTDIAG_WARMUP", 1)
	repetitions := fastDiagEnvInt(t, "FASTDIAG_REPETITIONS", 5)
	require.GreaterOrEqual(t, warmup, 0)
	require.GreaterOrEqual(t, repetitions, 1)
	require.LessOrEqual(t, warmup, 100)
	require.LessOrEqual(t, repetitions, 10, "repetitions must remain bounded for raw trace artifacts")

	params, residual, values, input := fastDiagP93Fixture(t)
	require.Equal(t, 55, bits.Len64(params.BootstrappingParameters.Q()[0]), "p93-q55 must use q0=55")
	inputSHA := fastDiagInputSHA256(values)
	require.Equal(t, fastDiagAcceptedInputSHA256, inputSHA, "p93-q55 workload must retain the accepted Count-1 fingerprint")
	baselineEval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	fastdiag.ConfigureCSV("")
	baseline, err := baselineEval.Bootstrap(input.CopyNew())
	require.NoError(t, err)
	wantDecoded := fastDiagDecode(t, residual, baseline)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)

	for i := 0; i < warmup; i++ {
		require.NoError(t, fastdiag.ConfigureCSV(scopes))
		fastdiag.Reset()
		_, err = eval.Bootstrap(input.CopyNew())
		require.NoError(t, err)
	}

	artifact := fastDiagTraceArtifact{
		SchemaVersion: "fastdiag.trace.v1", Timestamp: time.Now().UTC(), Profile: "p93-q55",
		Trace: strings.Split(scopes, ","), Warmup: warmup, Repetitions: repetitions,
		InputSHA256: inputSHA, InputLength: len(values), LogN: params.BootstrappingParameters.LogN(),
		LogSlots: params.CoeffsToSlotsParameters.LogSlots, QChainBits: fastDiagPrimeBits(params.BootstrappingParameters.Q()),
		PBits: fastDiagPrimeBits(params.BootstrappingParameters.P()), Polynomial: params.Mod1ParametersLiteral.Mod1Degree,
		DoubleAngle: params.Mod1ParametersLiteral.DoubleAngle, GoVersion: runtime.Version(), OS: runtime.GOOS,
		Arch: runtime.GOARCH, NumCPU: runtime.NumCPU(), GOMAXPROCS: runtime.GOMAXPROCS(0),
		Runs: make([]fastDiagRun, 0, repetitions),
	}
	for i := 0; i < repetitions; i++ {
		require.NoError(t, fastdiag.ConfigureCSV(scopes))
		fastdiag.Reset()
		started := time.Now()
		got, err := eval.Bootstrap(input.CopyNew())
		elapsed := time.Since(started).Nanoseconds()
		require.NoError(t, err)
		fastDiagRequireMetadataEqual(t, baseline, got)
		gotDecoded := fastDiagDecode(t, residual, got)
		require.Len(t, gotDecoded, len(wantDecoded))
		for slot := range gotDecoded {
			require.InDelta(t, real(wantDecoded[slot]), real(gotDecoded[slot]), 1e-2, "real slot %d", slot)
			require.InDelta(t, imag(wantDecoded[slot]), imag(gotDecoded[slot]), 1e-2, "imag slot %d", slot)
		}
		events := fastdiag.Events()
		fastDiagValidateEvents(t, scopes, events)
		artifact.Runs = append(artifact.Runs, fastDiagRun{Index: i + 1, ElapsedNS: elapsed, NumericalMatch: true, Events: events})
	}
	if path := os.Getenv("FASTDIAG_OUTPUT"); path != "" {
		data, err := json.MarshalIndent(artifact, "", "  ")
		require.NoError(t, err)
		require.NoError(t, os.WriteFile(path, append(data, '\n'), 0o600))
	}
}

func fastDiagInputSHA256(values []complex128) string {
	hash := sha256.New()
	var bytes [16]byte
	for _, value := range values {
		binary.LittleEndian.PutUint64(bytes[:8], math.Float64bits(real(value)))
		binary.LittleEndian.PutUint64(bytes[8:], math.Float64bits(imag(value)))
		_, _ = hash.Write(bytes[:])
	}
	return hex.EncodeToString(hash.Sum(nil))
}

func fastDiagDecode(t testing.TB, params ckks.Parameters, ct *rlwe.Ciphertext) []complex128 {
	t.Helper()
	encoder := ckks.NewEncoder(params)
	pt := ckks.NewPlaintext(params, ct.Level())
	*pt.MetaData = *ct.MetaData
	pt.Value.Copy(ct.Value[0])
	pt.IsNTT, pt.IsMontgomery = ct.IsNTT, ct.IsMontgomery
	decoded := make([]complex128, params.MaxSlots())
	require.NoError(t, encoder.Decode(pt, decoded))
	return decoded
}

func fastDiagRequireMetadataEqual(t testing.TB, want, got *rlwe.Ciphertext) {
	t.Helper()
	require.Equal(t, want.Level(), got.Level())
	require.Equal(t, want.Degree(), got.Degree())
	require.True(t, want.Scale.Equal(got.Scale))
	require.Equal(t, want.IsNTT, got.IsNTT)
	require.Equal(t, want.IsMontgomery, got.IsMontgomery)
	require.Equal(t, want.LogDimensions, got.LogDimensions)
}

func fastDiagValidateEvents(t testing.TB, scopes string, events []fastdiag.Event) {
	t.Helper()
	selected := map[string]bool{}
	for _, scope := range strings.Split(scopes, ",") {
		selected[strings.TrimSpace(scope)] = true
	}
	all := selected["all"]
	if selected["stage"] || all {
		var names []string
		for _, event := range events {
			if event.Scope == fastdiag.Stage {
				names = append(names, event.Name)
			}
		}
		require.Equal(t, []string{"bootstrap", "pack_n1_to_n2", "scale_down", "mod_up_trace", "coeffs_to_slots", "evalmod_real", "evalmod_imag", "slots_to_coeffs", "unpack_n2_to_n1", "public_finalization"}, names)
	}
	if selected["power"] || all {
		got := map[int]bool{}
		for _, event := range events {
			if event.Scope == fastdiag.Power && event.Name == "power" && event.Power != nil {
				got[*event.Power] = true
				require.NotNil(t, event.SplitA, "power %d is missing split A", *event.Power)
				require.NotNil(t, event.SplitB, "power %d is missing split B", *event.Power)
			}
		}
		for _, power := range []int{2, 3, 4, 6, 8, 16} {
			require.True(t, got[power], "missing generated power event %d", power)
		}
	}
	if selected["rescale"] || all {
		children := map[uint64][]fastdiag.Event{}
		for _, event := range events {
			if event.Scope == fastdiag.Rescale && event.ParentSequence != 0 {
				children[event.ParentSequence] = append(children[event.ParentSequence], event)
			}
		}
		parents := 0
		for _, event := range events {
			if event.Scope == fastdiag.Rescale && event.Name == "rescale" {
				parents++
				passes := children[event.Sequence]
				require.Equal(t, 2, len(passes), "Rescale %d must have exactly preflight and materialization children", event.Sequence)
				passByName := map[string]fastdiag.Event{}
				for _, pass := range passes {
					passByName[pass.Name] = pass
				}
				require.Len(t, passByName, 2, "Rescale %d has duplicate or unknown pass children", event.Sequence)
				preflight, hasPreflight := passByName["preflight"]
				materialization, hasMaterialization := passByName["materialization"]
				require.True(t, hasPreflight, "Rescale %d missing preflight child", event.Sequence)
				require.True(t, hasMaterialization, "Rescale %d missing materialization child", event.Sequence)
				fastDiagValidateRescalePass(t, preflight, children, false)
				fastDiagValidateRescalePass(t, materialization, children, true)
			}
		}
		require.Greater(t, parents, 0)
	}
}

func fastDiagValidateRescalePass(t testing.TB, pass fastdiag.Event, children map[uint64][]fastdiag.Event, materializing bool) {
	t.Helper()
	require.NotNil(t, pass.Count, "Rescale %s must record its component count", pass.Name)
	expectedComponents := *pass.Count
	passChildren := children[pass.Sequence]
	prefixByComponent := map[string]int{}
	restoreByComponent := map[string]int{}
	loopsByComponent := map[string]fastdiag.Event{}
	for _, child := range passChildren {
		switch child.Name {
		case "prefix_to_coefficient":
			prefixByComponent[child.Component]++
		case "coefficient_loop":
			require.NotContains(t, loopsByComponent, child.Component, "duplicate coefficient loop for %s/%s", pass.Name, child.Component)
			loopsByComponent[child.Component] = child
		case "ntt_montgomery_restore":
			require.True(t, materializing, "preflight must not contain an NTT restore event")
			restoreByComponent[child.Component]++
		default:
			t.Fatalf("unexpected direct child %q under Rescale %s", child.Name, pass.Name)
		}
	}
	if materializing {
		require.Empty(t, prefixByComponent, "commit must not repeat prefix conversion")
		require.Empty(t, loopsByComponent, "commit must not repeat coefficient reconstruction")
		require.Len(t, restoreByComponent, expectedComponents, "Rescale %s NTT restore component count", pass.Name)
		for component, count := range restoreByComponent {
			require.Equal(t, 1, count, "duplicate NTT restore for %s/%s", pass.Name, component)
		}
	} else {
		require.Len(t, prefixByComponent, expectedComponents, "Rescale %s prefix conversion component count", pass.Name)
		require.Len(t, loopsByComponent, expectedComponents, "Rescale %s coefficient-loop component count", pass.Name)
		require.Empty(t, restoreByComponent)
	}
	for component, count := range prefixByComponent {
		require.Equal(t, 1, count, "duplicate prefix conversion for %s/%s", pass.Name, component)
		loop, ok := loopsByComponent[component]
		require.True(t, ok, "Rescale %s component %s is missing its coefficient loop", pass.Name, component)
		loopChildren := children[loop.Sequence]
		wantNames := map[string]int{"reconstruct_center_round_capacity": 1, "residue_materialization": 1}
		gotNames := map[string]int{}
		for _, child := range loopChildren {
			gotNames[child.Name]++
		}
		require.Equal(t, wantNames, gotNames, "Rescale %s component %s coefficient-loop children", pass.Name, component)
	}
}

func fastDiagEnvInt(t testing.TB, name string, fallback int) int {
	t.Helper()
	value := os.Getenv(name)
	if value == "" {
		return fallback
	}
	parsed, err := strconv.Atoi(value)
	require.NoError(t, err, "%s must be an integer", name)
	return parsed
}

func fastDiagPrimeBits(primes []uint64) []int {
	out := make([]int, len(primes))
	for i, prime := range primes {
		out[i] = bits.Len64(prime)
	}
	return out
}
