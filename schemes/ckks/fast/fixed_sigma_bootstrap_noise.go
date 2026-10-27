package fast

import (
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// FixedSigmaBootstrapNoiseEvent is bounded per-multiplication evidence for
// the explicit research hook. The two sides are independently sampled even
// when the caller passed aliased operands.
type FixedSigmaBootstrapNoiseEvent struct {
	Operation         uint64                                  `json:"operation"`
	Level             int                                     `json:"level"`
	Scale             string                                  `json:"scale"`
	Left              ckks.FixedSigmaCoefficientNoiseEvidence `json:"left"`
	Right             ckks.FixedSigmaCoefficientNoiseEvidence `json:"right"`
	MultiplyCompleted bool                                    `json:"multiply_completed"`
	Failure           string                                  `json:"failure,omitempty"`
}

// FixedSigmaBootstrapNoiseSnapshot records the frozen configuration and the
// exact ct×ct calls that reached the injection boundary.
type FixedSigmaBootstrapNoiseSnapshot struct {
	Enabled             bool                            `json:"enabled"`
	Sigma               float64                         `json:"sigma_coefficient,omitempty"`
	RootSeed            int64                           `json:"root_seed,omitempty"`
	CtCtMultiplyCalls   uint64                          `json:"ct_ct_multiply_calls"`
	InjectedPolynomials uint64                          `json:"injected_polynomials"`
	Events              []FixedSigmaBootstrapNoiseEvent `json:"events,omitempty"`
}

type fixedSigmaBootstrapNoiseRuntime struct {
	sigma       float64
	rootSeed    int64
	operation   uint64
	helper      *ckks.Evaluator
	events      []FixedSigmaBootstrapNoiseEvent
	polynomials uint64
}

// EnableFixedSigmaBootstrapNoiseDiagnostic opts this Fast evaluator into the
// research-only pre-ct×ct hook. It is disabled by default. The coefficient
// sigma is fixed for the evaluator lifetime and is never calibrated here.
func (eval *Evaluator) EnableFixedSigmaBootstrapNoiseDiagnostic(sigma float64, rootSeed int64) error {
	if eval == nil {
		return fmt.Errorf("fixed-sigma Bootstrap noise evaluator cannot be nil")
	}
	if math.IsNaN(sigma) || math.IsInf(sigma, 0) || sigma <= 0 {
		return fmt.Errorf("fixed-sigma Bootstrap noise requires finite sigma > 0")
	}
	if err := validateFixedSigmaBootstrapNoiseParameters(eval.Parameters); err != nil {
		return fmt.Errorf("fixed-sigma Bootstrap noise CKKS parameters: %w", err)
	}
	// The existing coefficient-noise helper lives in the parent ckks package.
	// Fast already depends on that package for Parameters, so this reuses the
	// shared Gaussian/RNS/NTT and centered-capacity implementation without
	// changing the production evaluator or importing Fast back into its parent.
	help := ckks.NewEvaluator(eval.Parameters, nil)
	eval.fixedSigmaBootstrapNoise = &fixedSigmaBootstrapNoiseRuntime{
		sigma: sigma, rootSeed: rootSeed, helper: help,
	}
	return nil
}

// validateFixedSigmaBootstrapNoiseParameters checks structural and NTT
// capabilities required by the generic hook. Experiment profile restrictions
// (including a task's LogN) belong in the Primary measurement layer.
func validateFixedSigmaBootstrapNoiseParameters(params ckks.Parameters) error {
	n := params.N()
	logN := params.LogN()
	ringQ := params.RingQ()
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

// DisableFixedSigmaBootstrapNoiseDiagnostic restores the exact default-off
// route. It does not modify ciphertexts or reset any other evaluator state.
func (eval *Evaluator) DisableFixedSigmaBootstrapNoiseDiagnostic() {
	if eval != nil {
		eval.fixedSigmaBootstrapNoise = nil
	}
}

// FixedSigmaBootstrapNoiseDiagnosticSnapshot returns a copy of the evidence
// accumulated by the explicitly enabled hook.
func (eval *Evaluator) FixedSigmaBootstrapNoiseDiagnosticSnapshot() FixedSigmaBootstrapNoiseSnapshot {
	if eval == nil || eval.fixedSigmaBootstrapNoise == nil {
		return FixedSigmaBootstrapNoiseSnapshot{}
	}
	runtime := eval.fixedSigmaBootstrapNoise
	snapshot := FixedSigmaBootstrapNoiseSnapshot{
		Enabled: true, Sigma: runtime.sigma, RootSeed: runtime.rootSeed,
		CtCtMultiplyCalls: runtime.operation, InjectedPolynomials: runtime.polynomials,
		Events: append([]FixedSigmaBootstrapNoiseEvent(nil), runtime.events...),
	}
	for index := range snapshot.Events {
		snapshot.Events[index].Left.Phases = append([]ckks.FixedSigmaCoefficientNoisePhase(nil), runtime.events[index].Left.Phases...)
		snapshot.Events[index].Right.Phases = append([]ckks.FixedSigmaCoefficientNoisePhase(nil), runtime.events[index].Right.Phases...)
	}
	return snapshot
}

func (eval *Evaluator) applyFixedSigmaBootstrapNoise(
	op0 *rlwe.Ciphertext,
	op1 *rlwe.Element[ring.Poly],
	opOut *rlwe.Ciphertext,
	rows int,
	relin bool,
) error {
	runtime := eval.fixedSigmaBootstrapNoise
	runtime.operation++
	operation := runtime.operation
	left := op0.CopyNew()
	right := &rlwe.Ciphertext{Element: *op1.CopyNew()}
	leftSeed := fixedSigmaBootstrapNoiseSeed(runtime.rootSeed, operation, "left")
	rightSeed := fixedSigmaBootstrapNoiseSeed(runtime.rootSeed, operation, "right")
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
	err = eval.mulCore.ApplyRows(eval.Parameters.RingQ(), left, right.El(), opOut, rows, relin)
	if err != nil {
		event.Failure = fmt.Sprintf("Fast Mul ApplyRows: %v", err)
	} else {
		event.MultiplyCompleted = true
	}
	runtime.events = append(runtime.events, event)
	return err
}

func fixedSigmaBootstrapNoiseSeed(root int64, operation uint64, side string) int64 {
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
