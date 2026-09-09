package bootstrapping

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var _ Bootstrapper = (*FastEvaluator)(nil)

// ensureFastBootstrapCircuit lazily builds the Fast-only circuit data. The
// matrix and Mod1 planning is shared with Standard initialization, while all
// execution remains on the explicit Fast evaluators.
func (eval *FastEvaluator) ensureFastBootstrapCircuit() error {
	if eval == nil {
		return errors.New("Fast Bootstrap evaluator cannot be nil")
	}
	if eval.fastBootstrapInitialized {
		return eval.fastBootstrapErr
	}
	eval.fastBootstrapInitialized = true

	params := eval.Parameters
	if err := validateFastBootstrapParameters(params); err != nil {
		eval.fastBootstrapErr = err
		return err
	}
	adjusted, mod1Params, c2s, s2c, err := buildBootstrapCircuitData(params)
	if err != nil {
		eval.fastBootstrapErr = fmt.Errorf("cannot initialize Fast Bootstrap circuit: %w", err)
		return eval.fastBootstrapErr
	}
	eval.Parameters = adjusted
	eval.Mod1Parameters = mod1Params
	eval.C2SDFTMatrix = c2s
	eval.S2CDFTMatrix = s2c
	eval.DFTEvaluator = dftFastEvaluator(adjusted.BootstrappingParameters)
	eval.Mod1Evaluator = mod1.NewFastEvaluator(eval.FastCKKS, eval.PolynomialEvaluator, mod1Params)
	return nil
}

func dftFastEvaluator(params ckks.Parameters) *dft.FastEvaluator {
	return dft.NewFastEvaluator(params)
}

func validateFastBootstrapParameters(params Parameters) error {
	residual := params.ResidualParameters
	bootstrap := params.BootstrappingParameters
	if residual.RingType() != ring.Standard || bootstrap.RingType() != ring.Standard {
		return errors.New("Fast Bootstrap requires Standard rings")
	}
	if bootstrap.MaxLevel() < 1 {
		return errors.New("Fast Bootstrap requires at least q0 and q1 in BootstrappingParameters")
	}
	if residual.MaxLevel() > 1 {
		return fmt.Errorf("Fast Bootstrap Stage-A public output requires Residual MaxLevel <= 1, got %d", residual.MaxLevel())
	}
	if bootstrap.N() != residual.N() && bootstrap.N() != 2*residual.N() {
		return fmt.Errorf("Fast Bootstrap supports N2=N1 or N2=2*N1, got N1=%d N2=%d", residual.N(), bootstrap.N())
	}
	if params.CircuitOrder != ModUpThenEncode {
		return errors.New("Fast Bootstrap supports CircuitOrder ModUpThenEncode only")
	}
	if params.IterationsParameters != nil {
		return errors.New("Fast Bootstrap does not support iterative bootstrap parameters")
	}
	if params.EphemeralSecretWeight != 0 {
		return errors.New("Fast Bootstrap does not support ephemeral secret switching")
	}
	if params.Mod1ParametersLiteral.Mod1Type != mod1.CosDiscrete {
		return errors.New("Fast Bootstrap supports Mod1Type CosDiscrete only")
	}
	if params.Mod1ParametersLiteral.Mod1InvDegree != 0 {
		return errors.New("Fast Bootstrap does not support a Mod1 inverse polynomial")
	}
	if residual.RingQ() == nil || bootstrap.RingQ() == nil {
		return errors.New("Fast Bootstrap requires initialized residual and bootstrap parameters")
	}
	return nil
}

// validateFastBootstrapPublicInputs enforces the deliberately bounded public
// contract before any packing or in-place Fast operation can run.
func (eval *FastEvaluator) validateFastBootstrapPublicInputs(cts []rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast Bootstrap evaluator cannot be nil")
	}
	if len(cts) == 0 {
		return errors.New("Fast Bootstrap requires a non-empty ciphertext slice")
	}
	if err := validateFastBootstrapParameters(eval.Parameters); err != nil {
		return err
	}
	params := eval.Parameters.ResidualParameters
	var level, logSlots int
	for i := range cts {
		ct := &cts[i]
		if ct.MetaData == nil {
			return fmt.Errorf("ciphertext %d metadata cannot be nil", i)
		}
		if len(ct.Value) != 2 {
			return fmt.Errorf("ciphertext %d must have degree-one storage", i)
		}
		currentLevel := ct.Level()
		currentLogSlots := ct.LogSlots()
		if i == 0 {
			level, logSlots = currentLevel, currentLogSlots
			if logSlots < 0 || logSlots > eval.Parameters.LogMaxSlots() {
				return fmt.Errorf("ciphertext LogSlots %d is outside the supported packing dimensions", logSlots)
			}
		}
		if ct.N() != params.N() {
			return fmt.Errorf("ciphertext %d ring degree %d does not match ResidualParameters.N()=%d", i, ct.N(), params.N())
		}
		if ct.Degree() != 1 {
			return fmt.Errorf("ciphertext %d has degree %d, Fast Bootstrap requires degree 1", i, ct.Degree())
		}
		if currentLevel < 0 || currentLevel > params.MaxLevel() {
			return fmt.Errorf("ciphertext %d level %d is outside ResidualParameters", i, currentLevel)
		}
		if !ct.IsNTT {
			return fmt.Errorf("ciphertext %d is not in the NTT domain", i)
		}
		if ct.IsMontgomery {
			return fmt.Errorf("ciphertext %d is Montgomery at the public Fast Bootstrap boundary", i)
		}
		if ct.Scale.Cmp(rlwe.NewScale(0)) != 1 {
			return fmt.Errorf("ciphertext %d scale must be positive", i)
		}
		if len(ct.Value[0].Coeffs) < maintainedLimbs(currentLevel) || len(ct.Value[1].Coeffs) < maintainedLimbs(currentLevel) {
			return fmt.Errorf("ciphertext %d has invalid maintained storage", i)
		}
		if i > 0 {
			if currentLevel != level {
				return fmt.Errorf("ciphertext %d level %d does not match level %d", i, currentLevel, level)
			}
			if currentLogSlots != logSlots {
				return fmt.Errorf("ciphertext %d LogSlots %d does not match LogSlots %d", i, currentLogSlots, logSlots)
			}
			if ct.IsBatched != cts[0].IsBatched || ct.IsBitReversed != cts[0].IsBitReversed {
				return fmt.Errorf("ciphertext %d plaintext representation does not match ciphertext 0", i)
			}
		}
	}
	return nil
}

func (eval *FastEvaluator) bootstrapCore(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext, errScale *rlwe.Scale, err error) {
	if err = eval.ensureFastBootstrapCircuit(); err != nil {
		return nil, nil, err
	}
	if ctOut, errScale, err = eval.ScaleDown(ctIn); err != nil {
		return nil, nil, err
	}
	if ctOut, err = eval.ModUp(ctOut); err != nil {
		return nil, nil, err
	}
	var ctReal, ctImag *rlwe.Ciphertext
	if ctReal, ctImag, err = eval.DFTEvaluator.CoeffsToSlotsNew(ctOut, eval.C2SDFTMatrix); err != nil {
		return nil, nil, err
	}
	if ctReal, err = eval.EvalMod(ctReal); err != nil {
		return nil, nil, err
	}
	if ctImag != nil {
		if ctImag, err = eval.EvalMod(ctImag); err != nil {
			return nil, nil, err
		}
	}
	if ctOut, err = eval.DFTEvaluator.SlotsToCoeffsNew(ctReal, ctImag, eval.S2CDFTMatrix); err != nil {
		return nil, nil, err
	}
	return ctOut, errScale, nil
}

// Bootstrap applies the bounded Stage-A Fast bootstrap to one ciphertext.
func (eval *FastEvaluator) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if ct == nil {
		return nil, errors.New("Fast Bootstrap ciphertext cannot be nil")
	}
	cts, err := eval.BootstrapMany([]rlwe.Ciphertext{*ct})
	if err != nil {
		return nil, err
	}
	return &cts[0], nil
}

// BootstrapMany applies the public packing boundary, one Fast circuit per
// packed ciphertext, and the public Montgomery/scale restoration boundary.
func (eval *FastEvaluator) BootstrapMany(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error) {
	if err := eval.validateFastBootstrapPublicInputs(cts); err != nil {
		return nil, err
	}
	if err := eval.ensureFastBootstrapCircuit(); err != nil {
		return nil, err
	}
	packed, ctxtN1, ctxtN2, err := eval.PackAndSwitchN1ToN2(cts)
	if err != nil {
		return nil, fmt.Errorf("cannot Fast Bootstrap: %w", err)
	}
	for i := range packed {
		if packed[i].IsMontgomery {
			return nil, errors.New("Fast Bootstrap packing unexpectedly produced Montgomery public input")
		}
		if packed[i].Scale.Cmp(rlwe.NewScale(0)) != 1 {
			return nil, errors.New("Fast Bootstrap requires a positive ciphertext scale")
		}
		var coreOut *rlwe.Ciphertext
		if coreOut, _, err = eval.bootstrapCore(&packed[i]); err != nil {
			return nil, fmt.Errorf("cannot Fast Bootstrap circuit %d: %w", i, err)
		}
		packed[i] = *coreOut
	}
	if packed, err = eval.UnpackAndSwitchN2ToN1(packed, ctxtN1, ctxtN2); err != nil {
		return nil, fmt.Errorf("cannot Fast Bootstrap unpack: %w", err)
	}
	for i := range packed {
		if err = eval.finalizeFastPublicCiphertext(&packed[i]); err != nil {
			return nil, fmt.Errorf("cannot Fast Bootstrap finalize ciphertext %d: %w", i, err)
		}
	}
	return packed, nil
}

func (eval *FastEvaluator) finalizeFastPublicCiphertext(ct *rlwe.Ciphertext) error {
	if ct == nil || ct.MetaData == nil {
		return errors.New("Fast Bootstrap output and metadata cannot be nil")
	}
	params := eval.Parameters.ResidualParameters
	if ct.N() != params.N() || ct.Degree() != 1 || ct.Level() != params.MaxLevel() {
		return errors.New("Fast Bootstrap output does not match residual parameters")
	}
	if !ct.IsNTT || !ct.IsMontgomery {
		return errors.New("Fast Bootstrap output must be NTT-domain Montgomery before finalization")
	}
	ringQ := params.RingQ()
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < maintainedLimbs(ct.Level()); limb++ {
			ringQ.SubRings[limb].IMForm(ct.Value[d].Coeffs[limb], ct.Value[d].Coeffs[limb])
		}
	}
	ct.IsMontgomery = false
	ct.Scale = params.DefaultScale()
	return nil
}

func (eval *FastEvaluator) Depth() int {
	return eval.Parameters.BootstrappingParameters.MaxLevel() - eval.Parameters.ResidualParameters.MaxLevel()
}

func (eval *FastEvaluator) OutputLevel() int {
	return eval.Parameters.ResidualParameters.MaxLevel()
}

func (eval *FastEvaluator) MinimumInputLevel() int { return 0 }
