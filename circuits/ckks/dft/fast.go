package dft

import (
	"errors"
	"fmt"

	ltcommon "github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

// FastEvaluator executes the CKKS DFT factor sequence with the explicit
// q0/q1 Fast evaluator. It is intentionally separate from Evaluator: it does
// not implement the Standard LinearTransform, QP, or evaluation-key path.
type FastEvaluator struct {
	parameters ckks.Parameters
	eval       *fastckks.Evaluator
	scratch    fastDFTScratch
}

// fastDFTScratch keeps the ordinary ciphertext shape required by the public
// DFT surface, but allocates the two internal buffers only once. Fast DFT
// overwrites q0/q1 and never reads dormant limbs.
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
			copy: ckks.NewCiphertext(params, 1, params.MaxLevel()),
			tmp:  ckks.NewCiphertext(params, 1, params.MaxLevel()),
		},
	}
}

// FastEvaluator returns the underlying explicit Fast arithmetic evaluator.
func (eval *FastEvaluator) FastEvaluator() *fastckks.Evaluator { return eval.eval }

// CoeffsToSlotsNew applies Fast factorized CoeffsToSlots and returns its outputs.
func (eval *FastEvaluator) CoeffsToSlotsNew(ctIn *rlwe.Ciphertext, matrices Matrix) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	if err = eval.validateInput(ctIn); err != nil {
		return nil, nil, err
	}
	ctReal = ckks.NewCiphertext(eval.parameters, 1, matrices.LevelQ)
	setFastOutputDomain(ctReal, ctIn)
	if matrices.LogSlots == eval.parameters.LogMaxSlots() {
		ctImag = ckks.NewCiphertext(eval.parameters, 1, matrices.LevelQ)
		setFastOutputDomain(ctImag, ctIn)
	}
	err = eval.CoeffsToSlots(ctIn, matrices, ctReal, ctImag)
	return
}

// CoeffsToSlots applies Fast factorized CoeffsToSlots to the provided outputs.
func (eval *FastEvaluator) CoeffsToSlots(ctIn *rlwe.Ciphertext, matrices Matrix, ctReal, ctImag *rlwe.Ciphertext) error {
	if err := eval.validateInput(ctIn); err != nil {
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
		zV, err := eval.copyActive(ctIn)
		if err != nil {
			return err
		}
		if err = eval.dft(zV, matrices, zV); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots DFT: %w", err)
		}
		ctReal.Resize(1, zV.Level())
		if err = eval.eval.Conjugate(zV, ctReal); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots Conjugate: %w", err)
		}

		var tmp *rlwe.Ciphertext
		if ctImag != nil {
			ctImag.Resize(1, zV.Level())
			tmp = ctImag
		} else {
			tmp = eval.scratch.tmp
			tmp.Resize(1, zV.Level())
			setFastOutputDomain(tmp, zV)
		}
		if err = eval.eval.Sub(zV, ctReal, tmp); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots imaginary subtraction: %w", err)
		}
		if err = eval.eval.Mul(tmp, -1i, tmp); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots imaginary rotation: %w", err)
		}
		if err = eval.eval.Add(ctReal, zV, ctReal); err != nil {
			return fmt.Errorf("cannot Fast CoeffsToSlots real addition: %w", err)
		}
		if matrices.Format == RepackImagAsReal && matrices.LogSlots < eval.parameters.LogMaxSlots() {
			if err = eval.eval.Rotate(tmp, tmp, 1<<ctIn.LogDimensions.Cols); err != nil {
				return fmt.Errorf("cannot Fast CoeffsToSlots repack rotation: %w", err)
			}
			if err = eval.eval.Add(ctReal, tmp, ctReal); err != nil {
				return fmt.Errorf("cannot Fast CoeffsToSlots repack addition: %w", err)
			}
		}
		return nil
	}

	return eval.dft(ctIn, matrices, ctReal)
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
	opOut := ckks.NewCiphertext(eval.parameters, 1, matrices.LevelQ)
	setFastOutputDomain(opOut, ctReal)
	return opOut, eval.SlotsToCoeffs(ctReal, ctImag, matrices, opOut)
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
	setFastOutputDomain(opOut, ctReal)
	if ctImag != nil {
		if err := eval.eval.Mul(ctImag, 1i, opOut); err != nil {
			return fmt.Errorf("cannot Fast SlotsToCoeffs imaginary multiplication: %w", err)
		}
		if err := eval.eval.Add(opOut, ctReal, opOut); err != nil {
			return fmt.Errorf("cannot Fast SlotsToCoeffs real addition: %w", err)
		}
		return eval.dft(opOut, matrices, opOut)
	}
	return eval.dft(ctReal, matrices, opOut)
}

// dft evaluates each factor group and uses Fast Rescale once per group, which
// matches the Standard DFT factorization's logical level/scale progression.
func (eval *FastEvaluator) dft(ctIn *rlwe.Ciphertext, matrices Matrix, opOut *rlwe.Ciphertext) error {
	if len(matrices.Matrices) == 0 || len(matrices.Levels) == 0 {
		return errors.New("Fast DFT requires factor matrices")
	}
	inputDimensions := ctIn.LogDimensions
	matrixIdx := 0
	for _, factors := range matrices.Levels {
		for range factors {
			if matrixIdx >= len(matrices.Matrices) {
				return errors.New("Fast DFT factorization is shorter than its level schedule")
			}
			if err := eval.eval.LinearTransform(ctIn, ltcommon.LinearTransformation(matrices.Matrices[matrixIdx]), opOut); err != nil {
				return fmt.Errorf("Fast DFT factor %d: %w", matrixIdx, err)
			}
			matrixIdx++
			ctIn = opOut
		}
		if err := eval.eval.Rescale(opOut, opOut); err != nil {
			return fmt.Errorf("Fast DFT group rescale: %w", err)
		}
	}
	if matrixIdx != len(matrices.Matrices) {
		return errors.New("Fast DFT factorization has unused matrices")
	}
	opOut.LogDimensions = inputDimensions
	return nil
}

func (eval *FastEvaluator) copyActive(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if err := eval.validateInput(ct); err != nil {
		return nil, err
	}
	out := eval.scratch.copy
	out.Resize(1, ct.Level())
	*out.MetaData = *ct.MetaData
	for d := 0; d <= 1; d++ {
		copyQ01ForDFT(ct.Value[d], out.Value[d])
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

func copyQ01ForDFT(src, dst ring.Poly) {
	copy(dst.Coeffs[0], src.Coeffs[0])
	copy(dst.Coeffs[1], src.Coeffs[1])
}
