package ckks

import (
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// FixedSigmaBootstrapNoiseEvent records one public Fast CKKS ciphertext pair
// that reached the explicit pre-multiply noise hook.
type FixedSigmaBootstrapNoiseEvent struct {
	Operation         uint64                             `json:"operation"`
	Level             int                                `json:"level"`
	Scale             string                             `json:"scale"`
	Left              FixedSigmaCoefficientNoiseEvidence `json:"left"`
	Right             FixedSigmaCoefficientNoiseEvidence `json:"right"`
	MultiplyCompleted bool                               `json:"multiply_completed"`
	Failure           string                             `json:"failure,omitempty"`
}

// FixedSigmaBootstrapNoiseSnapshot is the bounded evidence accumulated by
// one explicitly configured CKKS evaluator. The hook is disabled by default.
type FixedSigmaBootstrapNoiseSnapshot struct {
	Enabled             bool                            `json:"enabled"`
	Sigma               float64                         `json:"sigma_coefficient,omitempty"`
	RootSeed            int64                           `json:"root_seed,omitempty"`
	CtCtMultiplyCalls   uint64                          `json:"ct_ct_multiply_calls"`
	InjectedPolynomials uint64                          `json:"injected_polynomials"`
	Events              []FixedSigmaBootstrapNoiseEvent `json:"events,omitempty"`
}

// FixedSigmaBootstrapNoiseRuntime contains the shared sampler/evidence path
// used by both the public Fast CKKS evaluator and the explicit Fast evaluator.
// Evaluators retain a pointer so their value-receiver operation methods share
// one deterministic counter and event log.
type FixedSigmaBootstrapNoiseRuntime struct {
	sigma       float64
	rootSeed    int64
	operation   uint64
	helper      *Evaluator
	events      []FixedSigmaBootstrapNoiseEvent
	polynomials uint64
}

// NewFixedSigmaBootstrapNoiseRuntime creates a disabled-by-default hook
// runtime for structurally valid Fast zero-secret CKKS parameters.
func NewFixedSigmaBootstrapNoiseRuntime(params Parameters, sigma float64, rootSeed int64) (*FixedSigmaBootstrapNoiseRuntime, error) {
	if math.IsNaN(sigma) || math.IsInf(sigma, 0) || sigma <= 0 {
		return nil, fmt.Errorf("fixed-sigma Bootstrap noise requires finite sigma > 0")
	}
	if err := validateFixedSigmaBootstrapNoiseParameters(params); err != nil {
		return nil, fmt.Errorf("fixed-sigma Bootstrap noise CKKS parameters: %w", err)
	}
	helper := NewEvaluator(params, nil)
	if !helper.usesFastCKKSZeroSecretSimulation() {
		return nil, fmt.Errorf("fixed-sigma Bootstrap noise is available only for Fast zero-secret CKKS parameters")
	}
	return &FixedSigmaBootstrapNoiseRuntime{sigma: sigma, rootSeed: rootSeed, helper: helper}, nil
}

func validateFixedSigmaBootstrapNoiseParameters(params Parameters) error {
	n, logN, ringQ := params.N(), params.LogN(), params.RingQ()
	if n < 2 || logN < 1 || ringQ == nil || len(ringQ.SubRings) == 0 || ringQ.SubRings[0] == nil {
		return fmt.Errorf("parameters do not contain an initialized CKKS ring")
	}
	if ringQ.N() != n || ringQ.LogN() != logN {
		return fmt.Errorf("CKKS N/LogN do not match RingQ (N=%d LogN=%d RingQ.N=%d RingQ.LogN=%d)", n, logN, ringQ.N(), ringQ.LogN())
	}
	if params.QCount() != len(ringQ.SubRings) || ringQ.Level() != params.MaxLevel() {
		return fmt.Errorf("CKKS Q chain does not match RingQ backing (QCount=%d subrings=%d maxLevel=%d ringLevel=%d)", params.QCount(), len(ringQ.SubRings), params.MaxLevel(), ringQ.Level())
	}
	if !params.NTTFlag() {
		return fmt.Errorf("CKKS parameters do not enable NTT representation")
	}
	for row, subring := range ringQ.SubRings {
		if subring == nil || subring.N != n {
			return fmt.Errorf("RingQ subring %d is missing or has an incompatible degree", row)
		}
	}
	return nil
}

// FixedSigmaBootstrapNoiseSeed derives the existing deterministic per-side
// seed from a root seed, operation index and operand side.
func FixedSigmaBootstrapNoiseSeed(root int64, operation uint64, side string) int64 {
	var input [17]byte
	binary.LittleEndian.PutUint64(input[0:8], uint64(root))
	binary.LittleEndian.PutUint64(input[8:16], operation)
	if side == "left" {
		input[16] = 0
	} else {
		input[16] = 1
	}
	digest := sha256.Sum256(input[:])
	return int64(binary.LittleEndian.Uint64(digest[:8]))
}

// Apply injects independent noise into copied c0 operands, runs exactly one
// caller-supplied multiplication, and records either its result or failure.
// The caller must route only validated Fast zero-secret ciphertext pairs here.
func (runtime *FixedSigmaBootstrapNoiseRuntime) Apply(
	op0 *rlwe.Ciphertext,
	op1 *rlwe.Element[ring.Poly],
	multiplyFailure string,
	multiply func(left *rlwe.Ciphertext, right *rlwe.Element[ring.Poly]) error,
) error {
	if runtime == nil || runtime.helper == nil || op0 == nil || op1 == nil || multiply == nil {
		return fmt.Errorf("fixed-sigma Bootstrap noise hook received an uninitialized multiplication")
	}
	runtime.operation++
	operation := runtime.operation
	left := op0.CopyNew()
	right := &rlwe.Ciphertext{Element: *op1.CopyNew()}
	leftSeed := FixedSigmaBootstrapNoiseSeed(runtime.rootSeed, operation, "left")
	rightSeed := FixedSigmaBootstrapNoiseSeed(runtime.rootSeed, operation, "right")
	event := FixedSigmaBootstrapNoiseEvent{Operation: operation, Level: min(op0.Level(), op1.Level()), Scale: op0.Scale.Value.Text('g', -1)}

	leftEvidence, err := runtime.helper.InjectFixedSigmaBootstrapMulNoiseDiagnostic(left, runtime.sigma, leftSeed)
	event.Left = leftEvidence
	if err != nil {
		event.Failure = fmt.Sprintf("left injection: %v", err)
		runtime.events = append(runtime.events, event)
		return err
	}
	rightEvidence, err := runtime.helper.InjectFixedSigmaBootstrapMulNoiseDiagnostic(right, runtime.sigma, rightSeed)
	event.Right = rightEvidence
	if err != nil {
		event.Failure = fmt.Sprintf("right injection: %v", err)
		runtime.events = append(runtime.events, event)
		return err
	}
	runtime.polynomials += 2
	err = multiply(left, right.El())
	if err != nil {
		if multiplyFailure == "" {
			multiplyFailure = "ct×ct multiply"
		}
		event.Failure = fmt.Sprintf("%s: %v", multiplyFailure, err)
	} else {
		event.MultiplyCompleted = true
	}
	runtime.events = append(runtime.events, event)
	return err
}

// Snapshot returns a defensive copy of all per-operation evidence.
func (runtime *FixedSigmaBootstrapNoiseRuntime) Snapshot() FixedSigmaBootstrapNoiseSnapshot {
	if runtime == nil {
		return FixedSigmaBootstrapNoiseSnapshot{}
	}
	snapshot := FixedSigmaBootstrapNoiseSnapshot{
		Enabled: true, Sigma: runtime.sigma, RootSeed: runtime.rootSeed,
		CtCtMultiplyCalls: runtime.operation, InjectedPolynomials: runtime.polynomials,
		Events: append([]FixedSigmaBootstrapNoiseEvent(nil), runtime.events...),
	}
	for index := range snapshot.Events {
		snapshot.Events[index].Left.Phases = append([]FixedSigmaCoefficientNoisePhase(nil), runtime.events[index].Left.Phases...)
		snapshot.Events[index].Right.Phases = append([]FixedSigmaCoefficientNoisePhase(nil), runtime.events[index].Right.Phases...)
	}
	return snapshot
}

// EnableFixedSigmaBootstrapNoiseDiagnostic explicitly enables the optional
// Fast-only hook on the public CKKS evaluator used by the original BTS.
func (eval *Evaluator) EnableFixedSigmaBootstrapNoiseDiagnostic(sigma float64, rootSeed int64) error {
	if eval == nil || eval.Encoder == nil || eval.Evaluator == nil {
		return fmt.Errorf("fixed-sigma Bootstrap noise evaluator cannot be nil or uninitialized")
	}
	runtime, err := NewFixedSigmaBootstrapNoiseRuntime(*eval.GetParameters(), sigma, rootSeed)
	if err != nil {
		return err
	}
	eval.fixedSigmaBootstrapNoise = runtime
	return nil
}

// DisableFixedSigmaBootstrapNoiseDiagnostic restores the exact default-off
// route without mutating ciphertexts or other evaluator state.
func (eval *Evaluator) DisableFixedSigmaBootstrapNoiseDiagnostic() {
	if eval != nil {
		eval.fixedSigmaBootstrapNoise = nil
	}
}

// FixedSigmaBootstrapNoiseDiagnosticSnapshot returns this evaluator's hook
// evidence, or an empty snapshot when the hook is disabled.
func (eval *Evaluator) FixedSigmaBootstrapNoiseDiagnosticSnapshot() FixedSigmaBootstrapNoiseSnapshot {
	if eval == nil {
		return FixedSigmaBootstrapNoiseSnapshot{}
	}
	return eval.fixedSigmaBootstrapNoise.Snapshot()
}
