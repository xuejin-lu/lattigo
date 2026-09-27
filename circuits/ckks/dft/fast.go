package dft

import (
	"errors"
	"fmt"
	"math/big"

	ltcommon "github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

// FastEvaluator executes the CKKS DFT factor sequence with the explicit
// Q-prefix Fast evaluator. It is intentionally separate from Evaluator: it does
// not implement the Standard LinearTransform, QP, or evaluation-key path.
type FastEvaluator struct {
	parameters ckks.Parameters
	eval       *fastckks.Evaluator
	scratch    fastDFTScratch
}

// fastDFTScratch keeps the ordinary ciphertext shape required by the public
// DFT surface, but allocates the two internal buffers only once. Fast DFT
// processes only the explicitly selected authoritative prefix.
type fastDFTScratch struct {
	copy *rlwe.Ciphertext
	tmp  *rlwe.Ciphertext
}

// NewFastEvaluator creates a bounded Fast DFT evaluator for Standard-ring
// NTT/Montgomery ciphertexts.
func NewFastEvaluator(params ckks.Parameters) *FastEvaluator {
	return &FastEvaluator{
		parameters: params,
		eval:       fastckks.NewEvaluator(params),
		scratch: fastDFTScratch{
			copy: fastckks.NewCiphertext(params, 1, params.MaxLevel()),
			tmp:  fastckks.NewCiphertext(params, 1, params.MaxLevel()),
		},
	}
}

// FastEvaluator returns the underlying explicit Fast arithmetic evaluator.
func (eval *FastEvaluator) FastEvaluator() *fastckks.Evaluator { return eval.eval }

// CoeffsToSlotsNew applies Fast factorized CoeffsToSlots and returns its outputs.
func (eval *FastEvaluator) CoeffsToSlotsNew(ctIn *rlwe.Ciphertext, matrices Matrix) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	return eval.CoeffsToSlotsNewWithRestorePlan(ctIn, matrices, nil)
}

// CoeffsToSlotsNewWithRestorePlan applies Fast factorized CoeffsToSlots and
// restores the maintained integer domain after the specified factor groups.
// A nil plan selects the ordinary Fast path.
func (eval *FastEvaluator) CoeffsToSlotsNewWithRestorePlan(ctIn *rlwe.Ciphertext, matrices Matrix, restorePlan []int) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	if err = eval.validateInput(ctIn); err != nil {
		return nil, nil, err
	}
	if err = validateRestorePlan(matrices, restorePlan); err != nil {
		return nil, nil, err
	}
	ctReal = fastckks.NewCiphertext(eval.parameters, 1, matrices.LevelQ)
	setFastOutputDomain(ctReal, ctIn)
	if matrices.LogSlots == eval.parameters.LogMaxSlots() {
		ctImag = fastckks.NewCiphertext(eval.parameters, 1, matrices.LevelQ)
		setFastOutputDomain(ctImag, ctIn)
	}
	err = eval.coeffsToSlots(ctIn, matrices, ctReal, ctImag, restorePlan)
	return
}

// CoeffsToSlots applies Fast factorized CoeffsToSlots to the provided outputs.
func (eval *FastEvaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext, matrices Matrix, ctReal, ctImag *rlwe.Ciphertext) error {
	return eval.coeffsToSlots(ctIn, matrices, ctReal, ctImag, nil)
}

// CoeffsToSlotsWithRestorePlan applies the bounded Fast C2S restore plan.
// Callers must provide one non-negative exponent per factor group.
func (eval *FastEvaluator) CoeffsToSlotsWithRestorePlan(ctIn *rlwe.Ciphertext, matrices Matrix, ctReal, ctImag *rlwe.Ciphertext, restorePlan []int) error {
	return eval.coeffsToSlots(ctIn, matrices, ctReal, ctImag, restorePlan)
}

func (eval *FastEvaluator) coeffsToSlots(ctIn *rlwe.Ciphertext, matrices Matrix, ctReal, ctImag *rlwe.Ciphertext, restorePlan []int) error {
	if err := eval.validateInput(ctIn); err != nil {
		return err
	}
	if err := validateRestorePlan(matrices, restorePlan); err != nil {
		return err
	}
	rows, err := fastckks.QPrefixWidth(ctIn.Level())
	if err != nil {
		return err
	}
	if ctReal == nil || matrices.Format == SplitRealAndImag && ctImag == nil && matrices.LogSlots == eval.parameters.LogMaxSlots() {
		return errors.New("Fast CoeffsToSlots requires the appropriate output ciphertexts")
	}
	setFastOutputDomain(ctReal, ctIn)
	if ctImag != nil {
		setFastOutputDomain(ctImag, ctIn)
	}

	if matrices.Format == RepackImagAsReal || matrices.Format == SplitRealAndImag {
		zV, err := eval.copyActive(ctIn, rows)
		if err != nil {
			return err
		}
		rows, err = eval.dftQPrefixRows(zV, matrices, zV, restorePlan, rows)
		if err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots DFT: %w", err)
		}
		fastckks.Resize(ctReal, 1, zV.Level(), eval.parameters.N())
		if err = eval.eval.ConjugateQPrefixRows(zV, ctReal, rows); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots Conjugate: %w", err)
		}

		var tmp *rlwe.Ciphertext
		if ctImag != nil {
			fastckks.Resize(ctImag, 1, zV.Level(), eval.parameters.N())
			tmp = ctImag
		} else {
			tmp = eval.scratch.tmp
			fastckks.Resize(tmp, 1, zV.Level(), eval.parameters.N())
			setFastOutputDomain(tmp, zV)
		}
		if err = eval.eval.SubQPrefixRows(zV, ctReal, tmp, rows); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots imaginary subtraction: %w", err)
		}
		if err = eval.eval.MulQPrefixRows(tmp, -1i, rows, tmp); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots imaginary rotation: %w", err)
		}
		if err = eval.eval.AddQPrefixRows(ctReal, zV, ctReal, rows); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots real addition: %w", err)
		}
		if matrices.Format == RepackImagAsReal && matrices.LogSlots < eval.parameters.LogMaxSlots() {
			if err = eval.eval.RotateQPrefixRows(tmp, tmp, 1<<ctIn.LogDimensions.Cols, rows); err != nil {
				return fmt.Errorf("cannot Fast CoeffsToSlots repack rotation: %w", err)
			}
			if err = eval.eval.AddQPrefixRows(ctReal, tmp, ctReal, rows); err != nil {
				return fmt.Errorf("cannot Fast CoeffsToSlots repack addition: %w", err)
			}
		}
		return nil
	}

	_, err = eval.dftQPrefixRows(ctIn, matrices, ctReal, restorePlan, rows)
	return err
}

// SlotsToCoeffsNew applies Fast factorized SlotsToCoeffs and returns its output.
func (eval *FastEvaluator) SlotsToCoeffsNew(ctReal, ctImag *rlwe.Ciphertext, matrices Matrix) (*rlwe.Ciphertext, error) {
	if err := eval.validateInput(ctReal); err != nil {
		return nil, err
	}
	if ctImag != nil {
		if err := eval.validateInput(ctImag); err != nil {
			return nil, err
		}
	}
	opOut := fastckks.NewCiphertext(eval.parameters, 1, matrices.LevelQ)
	setFastOutputDomain(opOut, ctReal)
	rows := fastckks.MaintainedLimbCount(&eval.parameters, ctReal.Level())
	return opOut, eval.SlotsToCoeffsQPrefixRows(ctReal, ctImag, matrices, opOut, rows)
}

// SlotsToCoeffsNewQPrefixRows applies S2C with an explicitly supplied input
// authority. Bootstrap uses this boundary to keep the current EvalMod
// producer's authority separate from its ciphertext backing width.
func (eval *FastEvaluator) SlotsToCoeffsNewQPrefixRows(ctReal, ctImag *rlwe.Ciphertext, matrices Matrix, rows int) (*rlwe.Ciphertext, error) {
	if err := eval.validateInput(ctReal); err != nil {
		return nil, err
	}
	if ctImag != nil {
		if err := eval.validateInput(ctImag); err != nil {
			return nil, err
		}
	}
	opOut := fastckks.NewCiphertext(eval.parameters, 1, matrices.LevelQ)
	setFastOutputDomain(opOut, ctReal)
	return opOut, eval.SlotsToCoeffsQPrefixRows(ctReal, ctImag, matrices, opOut, rows)
}

// SlotsToCoeffs applies Fast factorized SlotsToCoeffs to opOut.
func (eval *FastEvaluator) SlotsToCoeffs(ctReal, ctImag *rlwe.Ciphertext, matrices Matrix, opOut *rlwe.Ciphertext) error {
	if err := eval.validateInput(ctReal); err != nil {
		return err
	}
	if ctImag != nil {
		if err := eval.validateInput(ctImag); err != nil {
			return err
		}
	}
	if opOut == nil {
		return errors.New("Fast SlotsToCoeffs output cannot be nil")
	}
	rows := fastckks.MaintainedLimbCount(&eval.parameters, ctReal.Level())
	return eval.SlotsToCoeffsQPrefixRows(ctReal, ctImag, matrices, opOut, rows)
}

// SlotsToCoeffsQPrefixRows applies S2C using exactly rows explicit input
// authority, without promoting any physically allocated higher row.
func (eval *FastEvaluator) SlotsToCoeffsQPrefixRows(ctReal, ctImag *rlwe.Ciphertext, matrices Matrix, opOut *rlwe.Ciphertext, rows int) error {
	if err := eval.validateInput(ctReal); err != nil {
		return err
	}
	if ctImag != nil {
		if err := eval.validateInput(ctImag); err != nil {
			return err
		}
	}
	if opOut == nil {
		return errors.New("Fast SlotsToCoeffs output cannot be nil")
	}
	return eval.slotsToCoeffsQPrefixRows(ctReal, ctImag, matrices, opOut, rows)
}

func (eval *FastEvaluator) slotsToCoeffsQPrefixRows(ctReal, ctImag *rlwe.Ciphertext, matrices Matrix, opOut *rlwe.Ciphertext, rows int) error {
	if err := eval.validateInput(ctReal); err != nil {
		return err
	}
	if ctImag != nil {
		if err := eval.validateInput(ctImag); err != nil {
			return err
		}
		if ctImag.Level() != ctReal.Level() || !ctImag.Scale.Equal(ctReal.Scale) {
			return errors.New("Fast SlotsToCoeffs requires equal real/imaginary Levels and Scales")
		}
	}
	if opOut == nil {
		return errors.New("Fast SlotsToCoeffs output cannot be nil")
	}
	width, err := fastckks.QPrefixWidth(ctReal.Level())
	if err != nil {
		return err
	}
	if rows < 1 || rows > width {
		return fmt.Errorf("Fast SlotsToCoeffs row count %d must be in [1,%d] at level %d", rows, width, ctReal.Level())
	}
	fastckks.Resize(opOut, 1, ctReal.Level(), eval.parameters.N())
	setFastOutputDomain(opOut, ctReal)
	if ctImag != nil {
		if err := eval.eval.MulQPrefixRows(ctImag, 1i, rows, opOut); err != nil {
			return fmt.Errorf("cannot Fast SlotsToCoeffs imaginary multiplication: %w", err)
		}
		if err := eval.eval.AddQPrefixRows(opOut, ctReal, opOut, rows); err != nil {
			return fmt.Errorf("cannot Fast SlotsToCoeffs real addition: %w", err)
		}
		_, err := eval.dftQPrefixRows(opOut, matrices, opOut, nil, rows)
		return err
	}
	_, err = eval.dftQPrefixRows(ctReal, matrices, opOut, nil, rows)
	return err
}

// dft evaluates each factor group and uses Fast Rescale once per group, which
// matches the Standard DFT factorization's logical level/scale progression.
func (eval *FastEvaluator) dft(ctIn *rlwe.Ciphertext, matrices Matrix, opOut *rlwe.Ciphertext) error {
	rows := fastckks.MaintainedLimbCount(&eval.parameters, ctIn.Level())
	_, err := eval.dftQPrefixRows(ctIn, matrices, opOut, nil, rows)
	return err
}

func (eval *FastEvaluator) dftWithRestorePlan(ctIn *rlwe.Ciphertext, matrices Matrix, opOut *rlwe.Ciphertext, restorePlan []int) error {
	rows := fastckks.MaintainedLimbCount(&eval.parameters, ctIn.Level())
	_, err := eval.dftQPrefixRows(ctIn, matrices, opOut, restorePlan, rows)
	return err
}

// dftQPrefixRows executes factor groups using an explicit input authority.
// After each Rescale, authority contracts to the policy width at the new
// logical Level and is never inferred from backing allocation.
func (eval *FastEvaluator) dftQPrefixRows(ctIn *rlwe.Ciphertext, matrices Matrix, opOut *rlwe.Ciphertext, restorePlan []int, initialRows int) (rows int, err error) {
	if len(matrices.Matrices) == 0 || len(matrices.Levels) == 0 {
		return 0, errors.New("Fast DFT requires factor matrices")
	}
	if err := validateRestorePlan(matrices, restorePlan); err != nil {
		return 0, err
	}
	if err := eval.validateInput(ctIn); err != nil {
		return 0, err
	}
	if opOut == nil {
		return 0, errors.New("Fast DFT output cannot be nil")
	}
	width, err := fastckks.QPrefixWidth(ctIn.Level())
	if err != nil {
		return 0, err
	}
	if initialRows < 1 || initialRows > width {
		return 0, fmt.Errorf("Fast DFT input row count %d must be in [1,%d] at level %d", initialRows, width, ctIn.Level())
	}
	rows = initialRows
	inputDimensions := ctIn.LogDimensions
	matrixIdx := 0
	for groupIdx, factors := range matrices.Levels {
		for range factors {
			if matrixIdx >= len(matrices.Matrices) {
				return rows, errors.New("Fast DFT factorization is shorter than its level schedule")
			}
			if err := eval.eval.LinearTransformQPrefixRows(ctIn, ltcommon.LinearTransformation(matrices.Matrices[matrixIdx]), rows, opOut); err != nil {
				return rows, fmt.Errorf("Fast DFT factor %d: %w", matrixIdx, err)
			}
			matrixIdx++
			ctIn = opOut
		}
		if err := eval.eval.RescaleQPrefixRows(opOut, rows, opOut); err != nil {
			return rows, fmt.Errorf("Fast DFT group rescale: %w", err)
		}
		targetRows, err := fastckks.QPrefixWidth(opOut.Level())
		if err != nil {
			return rows, err
		}
		rows = min(rows, targetRows)
		if restorePlan != nil && restorePlan[groupIdx] != 0 {
			k := restorePlan[groupIdx]
			scalar := new(big.Int).Lsh(big.NewInt(1), uint(k))
			if err := eval.eval.MulIntegerQPrefixRows(opOut, scalar, rows, opOut); err != nil {
				return rows, fmt.Errorf("Fast DFT group %d restore: %w", groupIdx, err)
			}
			opOut.Scale = opOut.Scale.Mul(rlwe.NewScale(scalar))
		}
	}
	if matrixIdx != len(matrices.Matrices) {
		return rows, errors.New("Fast DFT factorization has unused matrices")
	}
	opOut.LogDimensions = inputDimensions
	return rows, nil
}

func validateRestorePlan(matrices Matrix, restorePlan []int) error {
	if restorePlan == nil {
		return nil
	}
	if len(restorePlan) != len(matrices.Levels) {
		return fmt.Errorf("Fast DFT restore plan has %d groups, expected %d", len(restorePlan), len(matrices.Levels))
	}
	for i, k := range restorePlan {
		if k < 0 || k > 30 {
			return fmt.Errorf("Fast DFT restore exponent %d at group %d is unsupported", k, i)
		}
	}
	return nil
}

func (eval *FastEvaluator) copyActive(ct *rlwe.Ciphertext, rows int) (*rlwe.Ciphertext, error) {
	if err := eval.validateInput(ct); err != nil {
		return nil, err
	}
	width, err := fastckks.QPrefixWidth(ct.Level())
	if err != nil {
		return nil, err
	}
	if rows < 1 || rows > width {
		return nil, fmt.Errorf("Fast DFT copy row count %d must be in [1,%d] at level %d", rows, width, ct.Level())
	}
	for component := 0; component <= 1; component++ {
		if len(ct.Value[component].Coeffs) < rows {
			return nil, fmt.Errorf("Fast DFT input component %d has %d rows, requested %d", component, len(ct.Value[component].Coeffs), rows)
		}
		for row := 0; row < rows; row++ {
			if len(ct.Value[component].Coeffs[row]) != eval.parameters.N() {
				return nil, fmt.Errorf("Fast DFT input component %d q%d backing is invalid", component, row)
			}
		}
	}
	out := eval.scratch.copy
	fastckks.Resize(out, 1, ct.Level(), eval.parameters.N())
	*out.MetaData = *ct.MetaData
	for d := 0; d <= 1; d++ {
		copyPrefixRowsForDFT(rows, ct.Value[d], out.Value[d])
	}
	return out, nil
}

func (eval *FastEvaluator) validateInput(ct *rlwe.Ciphertext) error {
	if eval == nil || eval.eval == nil || ct == nil || ct.MetaData == nil {
		return errors.New("Fast DFT evaluator and ciphertext metadata cannot be nil")
	}
	if ct.N() != eval.parameters.N() || ct.Degree() != 1 || ct.Level() < 1 {
		return errors.New("Fast DFT requires a degree-one ciphertext with q0 and q1")
	}
	if eval.parameters.RingType() != ring.Standard {
		return errors.New("Fast DFT requires the Standard ring")
	}
	if !ct.IsNTT || !ct.IsMontgomery {
		return errors.New("Fast DFT requires NTT/Montgomery ciphertexts")
	}
	return nil
}

func setFastOutputDomain(out, reference *rlwe.Ciphertext) {
	out.IsNTT = reference.IsNTT
	out.IsMontgomery = reference.IsMontgomery
	out.IsBatched = reference.IsBatched
	out.IsBitReversed = reference.IsBitReversed
}

func copyPrefixRowsForDFT(rows int, src, dst ring.Poly) {
	for row := 0; row < rows; row++ {
		copy(dst.Coeffs[row], src.Coeffs[row])
	}
}
