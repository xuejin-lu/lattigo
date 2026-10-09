package bootstrapping

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/internal/fastdiag"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
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
	if c2s, eval.C2SRestorePlan, eval.C2SCompressionActive, err = prepareFastLogN13C2S(adjusted, c2s); err != nil {
		eval.fastBootstrapErr = fmt.Errorf("cannot prepare Fast LogN13 C2S compression: %w", err)
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
	// Preserve the Standard frontend's E=32 value. Fast Bootstrap simulates a
	// zero secret and does not execute the dense/sparse secret KeySwitch path.
	if params.EphemeralSecretWeight != 0 && params.EphemeralSecretWeight != 32 {
		return fmt.Errorf("Fast Bootstrap supports EphemeralSecretWeight 0 or 32, got %d", params.EphemeralSecretWeight)
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
		rows, err := fastckks.QPrefixWidth(currentLevel)
		if err != nil {
			return fmt.Errorf("ciphertext %d: %w", i, err)
		}
		for d := 0; d <= 1; d++ {
			if len(ct.Value[d].Coeffs) < rows {
				return fmt.Errorf("ciphertext %d component %d has insufficient public Q-prefix rows", i, d)
			}
			for row := 0; row < rows; row++ {
				if len(ct.Value[d].Coeffs[row]) != params.N() {
					return fmt.Errorf("ciphertext %d component %d q%d row has invalid public backing", i, d, row)
				}
			}
		}
		for row := 0; row < rows; row++ {
			for coefficient, value := range ct.Value[1].Coeffs[row] {
				if value != 0 {
					return fmt.Errorf("ciphertext %d c1 q%d coefficient %d is non-zero; Fast Bootstrap requires a zero-secret simulation input", i, row, coefficient)
				}
			}
		}
		if i > 0 {
			if currentLevel != level {
				return fmt.Errorf("ciphertext %d level %d does not match level %d", i, currentLevel, level)
			}
			if currentLogSlots != logSlots {
				return fmt.Errorf("ciphertext %d LogSlots %d does not match LogSlots %d", i, currentLogSlots, logSlots)
			}
			if ct.IsNTT != cts[0].IsNTT || ct.IsMontgomery != cts[0].IsMontgomery ||
				ct.IsBatched != cts[0].IsBatched || ct.IsBitReversed != cts[0].IsBitReversed {
				return fmt.Errorf("ciphertext %d plaintext representation does not match ciphertext 0", i)
			}
			if !ct.Scale.Equal(cts[0].Scale) || ct.LogDimensions != cts[0].LogDimensions {
				return fmt.Errorf("ciphertext %d scale or dimensions do not match ciphertext 0", i)
			}
		}
	}
	return nil
}

func fastDiagOutput(ct *rlwe.Ciphertext) fastdiag.Fields {
	if ct == nil {
		return fastdiag.Fields{}
	}
	rows, _ := fastckks.QPrefixWidth(ct.Level())
	return fastdiag.Output(ct, rows)
}

func (eval *FastEvaluator) bootstrapCore(ctIn *rlwe.Ciphertext, parent uint64) (ctOut *rlwe.Ciphertext, errScale *rlwe.Scale, err error) {
	if err = eval.ensureFastBootstrapCircuit(); err != nil {
		return nil, nil, err
	}
	var stageSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan = fastdiag.Begin(fastdiag.Stage, "scale_down", parent, fastdiag.Input(ctIn, -1))
	}
	if ctOut, errScale, err = eval.ScaleDown(ctIn); err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan.End(fastdiag.Fields{})
		}
		return nil, nil, err
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan.End(fastDiagOutput(ctOut))
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan = fastdiag.Begin(fastdiag.Stage, "mod_up_trace", parent, fastdiag.Input(ctOut, -1))
	}
	if ctOut, err = eval.ModUp(ctOut); err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan.End(fastdiag.Fields{})
		}
		return nil, nil, err
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan.End(fastDiagOutput(ctOut))
	}
	var ctReal, ctImag *rlwe.Ciphertext
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan = fastdiag.Begin(fastdiag.Stage, "coeffs_to_slots", parent, fastdiag.Input(ctOut, -1).Merge(fastdiag.RepetitionCount(1)))
	}
	if ctReal, ctImag, err = eval.DFTEvaluator.CoeffsToSlotsNewWithRestorePlan(ctOut, eval.C2SDFTMatrix, eval.C2SRestorePlan); err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan.End(fastdiag.Fields{})
		}
		return nil, nil, err
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan.End(fastDiagOutput(ctReal).Merge(fastdiag.ComponentName("real+imag")))
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan = fastdiag.Begin(fastdiag.Stage, "evalmod_real", parent, fastdiag.Input(ctReal, -1).Merge(fastdiag.ComponentName("real")))
	}
	if ctReal, err = eval.EvalMod(ctReal); err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan.End(fastdiag.Fields{})
		}
		return nil, nil, err
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan.End(fastDiagOutput(ctReal))
	}
	if ctImag != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan = fastdiag.Begin(fastdiag.Stage, "evalmod_imag", parent, fastdiag.Input(ctImag, -1).Merge(fastdiag.ComponentName("imag")))
		}
		if ctImag, err = eval.EvalMod(ctImag); err != nil {
			if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
				stageSpan.End(fastdiag.Fields{})
			}
			return nil, nil, err
		}
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan.End(fastDiagOutput(ctImag))
		}
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan = fastdiag.Begin(fastdiag.Stage, "slots_to_coeffs", parent, fastdiag.Input(ctReal, -1).Merge(fastdiag.ComponentName("real+imag")))
	}
	if ctOut, err = eval.SlotsToCoeffs(ctReal, ctImag); err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan.End(fastdiag.Fields{})
		}
		return nil, nil, err
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan.End(fastDiagOutput(ctOut))
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
	var bootstrapSpan fastdiag.Span
	var parent uint64
	var traceOutput *rlwe.Ciphertext
	var traceOutputRows = -1
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		var input *rlwe.Ciphertext
		if len(cts) != 0 {
			input = &cts[0]
		}
		bootstrapSpan = fastdiag.Begin(fastdiag.Stage, "bootstrap", 0, fastdiag.Input(input, -1))
		parent = bootstrapSpan.Sequence()
	}
	if err := eval.validateFastBootstrapPublicInputs(cts); err != nil {
		return nil, err
	}
	if err := eval.ensureFastBootstrapCircuit(); err != nil {
		return nil, err
	}
	inputRows, err := fastckks.QPrefixWidth(cts[0].Level())
	if err != nil {
		return nil, err
	}
	var stageSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		var input *rlwe.Ciphertext
		if len(cts) != 0 {
			input = &cts[0]
		}
		stageSpan = fastdiag.Begin(fastdiag.Stage, "pack_n1_to_n2", parent, fastdiag.Input(input, inputRows))
	}
	packed, ctxtN1, ctxtN2, err := eval.PackAndSwitchN1ToN2QPrefixRows(cts, inputRows)
	if err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan.End(fastdiag.Fields{})
		}
		return nil, fmt.Errorf("cannot Fast Bootstrap: %w", err)
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		var output *rlwe.Ciphertext
		if len(packed) != 0 {
			output = &packed[0]
		}
		stageSpan.End(fastDiagOutput(output).Merge(fastdiag.RepetitionCount(len(packed))))
	}
	for i := range packed {
		if packed[i].IsMontgomery {
			return nil, errors.New("Fast Bootstrap packing unexpectedly produced Montgomery public input")
		}
		if packed[i].Scale.Cmp(rlwe.NewScale(0)) != 1 {
			return nil, errors.New("Fast Bootstrap requires a positive ciphertext scale")
		}
		var coreOut *rlwe.Ciphertext
		if coreOut, _, err = eval.bootstrapCore(&packed[i], parent); err != nil {
			return nil, fmt.Errorf("cannot Fast Bootstrap circuit %d: %w", i, err)
		}
		packed[i] = *coreOut
	}
	outputRows, err := fastckks.QPrefixWidth(packed[0].Level())
	if err != nil {
		return nil, err
	}
	if outputRows > eval.Parameters.ResidualParameters.MaxLevel()+1 {
		return nil, fmt.Errorf("Fast Bootstrap core output Q-prefix width %d exceeds the public residual width", outputRows)
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan = fastdiag.Begin(fastdiag.Stage, "unpack_n2_to_n1", parent, fastdiag.Input(&packed[0], outputRows))
	}
	if packed, err = eval.UnpackAndSwitchN2ToN1QPrefixRows(packed, ctxtN1, ctxtN2, outputRows); err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
			stageSpan.End(fastdiag.Fields{})
		}
		return nil, fmt.Errorf("cannot Fast Bootstrap unpack: %w", err)
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		var output *rlwe.Ciphertext
		if len(packed) != 0 {
			output = &packed[0]
		}
		rows, _ := fastckks.QPrefixWidth(output.Level())
		stageSpan.End(fastDiagOutput(output).Merge(fastdiag.RepetitionCount(len(packed))).Merge(fastdiag.Fields{RowsOut: &rows}))
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		stageSpan = fastdiag.Begin(fastdiag.Stage, "public_finalization", parent, fastdiag.Input(&packed[0], outputRows))
	}
	for i := range packed {
		if err = eval.finalizeFastPublicCiphertext(&packed[i]); err != nil {
			if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
				stageSpan.End(fastdiag.Fields{})
			}
			return nil, fmt.Errorf("cannot Fast Bootstrap finalize ciphertext %d: %w", i, err)
		}
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		var output *rlwe.Ciphertext
		if len(packed) != 0 {
			output = &packed[0]
		}
		traceOutput = output
		if output != nil {
			traceOutputRows, _ = fastckks.QPrefixWidth(output.Level())
		}
		stageSpan.End(fastDiagOutput(output).Merge(fastdiag.RepetitionCount(len(packed))))
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Stage) {
		bootstrapSpan.End(fastdiag.Output(traceOutput, traceOutputRows))
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
	rows, err := fastckks.QPrefixWidth(ct.Level())
	if err != nil {
		return err
	}
	if rows > params.MaxLevel()+1 {
		return errors.New("Fast Bootstrap output authority exceeds the public residual Level")
	}
	for d := 0; d <= 1; d++ {
		if len(ct.Value[d].Coeffs) < rows {
			return fmt.Errorf("Fast Bootstrap output component %d has insufficient public rows", d)
		}
		for row := 0; row < rows; row++ {
			if len(ct.Value[d].Coeffs[row]) != params.N() {
				return fmt.Errorf("Fast Bootstrap output component %d q%d row has invalid public backing", d, row)
			}
		}
	}
	ringQ := params.RingQ()
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < rows; limb++ {
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
