package fast

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// FixedSigmaBootstrapNoiseEvent and FixedSigmaBootstrapNoiseSnapshot are
// aliases so the explicit Fast evaluator and public CKKS evaluator share the
// same hook evidence schema.
type FixedSigmaBootstrapNoiseEvent = ckks.FixedSigmaBootstrapNoiseEvent
type FixedSigmaBootstrapNoiseSnapshot = ckks.FixedSigmaBootstrapNoiseSnapshot
type fixedSigmaBootstrapNoiseRuntime = ckks.FixedSigmaBootstrapNoiseRuntime

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
	runtime, err := ckks.NewFixedSigmaBootstrapNoiseRuntime(eval.Parameters, sigma, rootSeed)
	if err != nil {
		return err
	}
	eval.fixedSigmaBootstrapNoise = runtime
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
	return eval.fixedSigmaBootstrapNoise.Snapshot()
}

func (eval *Evaluator) applyFixedSigmaBootstrapNoise(
	op0 *rlwe.Ciphertext,
	op1 *rlwe.Element[ring.Poly],
	opOut *rlwe.Ciphertext,
	rows int,
	relin bool,
) error {
	return eval.fixedSigmaBootstrapNoise.Apply(op0, op1, "Fast Mul ApplyRows", func(left *rlwe.Ciphertext, right *rlwe.Element[ring.Poly]) error {
		return eval.mulCore.ApplyRows(eval.Parameters.RingQ(), left, right, opOut, rows, relin)
	})
}
