package fastcore

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// ScalarCore is the import-neutral shared boundary for the existing Fast
// complex-scalar kernels. Both the ordinary CKKS adapter and explicit Fast
// evaluator use the same Q-prefix arithmetic implementation.
type ScalarCore interface {
	AddSubRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows int, scalar *bignum.Complex, sub bool) error
	MulRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows, levelsConsumed int, scalar *bignum.Complex) error
	MulAtScaleRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows int, scalar *bignum.Complex, scale rlwe.Scale) error
	MulThenAddRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows, levelsConsumed int, scalar *bignum.Complex) error
	MulThenAddEncodedRows(ringQ *ring.Ring, source []ring.Poly, opOut *rlwe.Ciphertext, rows int, values [MaxQPrefixWidth][2]uint64) error
}

// ScalarWorkspace shares the Fast complex-scalar arithmetic and keeps the
// bounded scratch needed when an in-place MulThenAdd must first align scales.
// Like the evaluator workspaces, it is intended for one execution stream.
type ScalarWorkspace struct {
	n       int
	maxRows int
	scratch []ring.Poly
}

func NewScalarWorkspace(ringQ *ring.Ring, maxRows int) *ScalarWorkspace {
	w := &ScalarWorkspace{maxRows: maxRows}
	if ringQ == nil || maxRows < 1 || maxRows > MaxQPrefixWidth {
		return w
	}
	w.n = ringQ.N()
	w.scratch = make([]ring.Poly, 3)
	for i := range w.scratch {
		w.scratch[i] = ring.NewPoly(w.n, maxRows-1)
	}
	return w
}

// ScalarNTTRows encodes a complex scalar into the two NTT lanes for exactly
// rows authoritative Q-prefix moduli. It is the common implementation used
// by both public and explicit Fast CKKS evaluators.
func ScalarNTTRows(ringQ *ring.Ring, scalar *bignum.Complex, scale *big.Float, montgomery bool, level, rows int) ([MaxQPrefixWidth][2]uint64, error) {
	var out [MaxQPrefixWidth][2]uint64
	if ringQ == nil || scalar == nil || scale == nil {
		return out, errors.New("Fast scalar ring, scalar and scale cannot be nil")
	}
	if err := ValidatePrefixRows(ringQ, level, rows); err != nil {
		return out, err
	}
	toInt := func(x *big.Float) *big.Int {
		z := new(big.Int)
		if x == nil {
			return z
		}
		v := new(big.Float).Mul(x, scale)
		if x.Sign() > 0 {
			v.Add(v, new(big.Float).SetFloat64(0.5))
		} else if x.Sign() < 0 {
			v.Sub(v, new(big.Float).SetFloat64(0.5))
		}
		v.Int(z)
		return z
	}
	real, imag := toInt(scalar[0]), toInt(scalar[1])
	for limb := 0; limb < rows; limb++ {
		subring := ringQ.SubRings[limb]
		r := new(big.Int).Mod(new(big.Int).Set(real), new(big.Int).SetUint64(subring.Modulus)).Uint64()
		i := new(big.Int).Mod(new(big.Int).Set(imag), new(big.Int).SetUint64(subring.Modulus)).Uint64()
		i = ring.MRed(i, subring.RootsForward[1], subring.Modulus, subring.MRedConstant)
		out[limb][0], out[limb][1] = ring.CRed(r+i, subring.Modulus), ring.CRed(r+subring.Modulus-i, subring.Modulus)
		if montgomery {
			out[limb][0] = ring.MForm(out[limb][0], subring.Modulus, subring.BRedConstant)
			out[limb][1] = ring.MForm(out[limb][1], subring.Modulus, subring.BRedConstant)
		}
	}
	return out, nil
}

func (w *ScalarWorkspace) AddSubRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows int, scalar *bignum.Complex, sub bool) error {
	level, err := w.validateUnary(ringQ, op0, opOut, rows)
	if err != nil {
		return err
	}
	if scalar == nil {
		return errors.New("Fast scalar Add/Sub scalar cannot be nil")
	}
	values, err := ScalarNTTRows(ringQ, scalar, &op0.Scale.Value, op0.IsMontgomery, level, rows)
	if err != nil {
		return err
	}
	w.resizeOutput(op0, opOut, level)
	if op0 != opOut {
		for component := range op0.Value {
			CopyPrefixRowsUnchecked(rows, op0.Value[component], opOut.Value[component])
		}
	}
	for limb := 0; limb < rows; limb++ {
		subring := ringQ.SubRings[limb]
		half := w.n >> 1
		if sub {
			subring.SubScalar(opOut.Value[0].Coeffs[limb][:half], values[limb][0], opOut.Value[0].Coeffs[limb][:half])
			subring.SubScalar(opOut.Value[0].Coeffs[limb][half:], values[limb][1], opOut.Value[0].Coeffs[limb][half:])
		} else {
			subring.AddScalar(opOut.Value[0].Coeffs[limb][:half], values[limb][0], opOut.Value[0].Coeffs[limb][:half])
			subring.AddScalar(opOut.Value[0].Coeffs[limb][half:], values[limb][1], opOut.Value[0].Coeffs[limb][half:])
		}
	}
	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale
	opOut.IsNTT, opOut.IsMontgomery, opOut.IsBatched = op0.IsNTT, op0.IsMontgomery, op0.IsBatched
	return nil
}

func (w *ScalarWorkspace) MulRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows, levelsConsumed int, scalar *bignum.Complex) error {
	level, err := w.validateUnary(ringQ, op0, opOut, rows)
	if err != nil {
		return err
	}
	if scalar == nil {
		return errors.New("Fast scalar Mul scalar cannot be nil")
	}
	scale := rlwe.NewScale(1)
	if !scalar.IsInt() {
		scale, err = CoefficientScale(ringQ, level, levelsConsumed)
		if err != nil {
			return err
		}
	}
	return w.MulAtScaleRows(ringQ, op0, opOut, rows, scalar, scale)
}

func (w *ScalarWorkspace) MulAtScaleRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows int, scalar *bignum.Complex, scale rlwe.Scale) error {
	level, err := w.validateUnary(ringQ, op0, opOut, rows)
	if err != nil {
		return err
	}
	if scalar == nil {
		return errors.New("Fast scalar Mul scalar cannot be nil")
	}
	values, err := ScalarNTTRows(ringQ, scalar, &scale.Value, false, level, rows)
	if err != nil {
		return err
	}
	w.resizeOutput(op0, opOut, level)
	for component := range op0.Value {
		for limb := 0; limb < rows; limb++ {
			subring := ringQ.SubRings[limb]
			half := w.n >> 1
			scalar0 := ring.MForm(values[limb][0], subring.Modulus, subring.BRedConstant)
			scalar1 := ring.MForm(values[limb][1], subring.Modulus, subring.BRedConstant)
			subring.MulScalarMontgomery(op0.Value[component].Coeffs[limb][:half], scalar0, opOut.Value[component].Coeffs[limb][:half])
			subring.MulScalarMontgomery(op0.Value[component].Coeffs[limb][half:], scalar1, opOut.Value[component].Coeffs[limb][half:])
		}
	}
	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale.Mul(scale)
	opOut.IsNTT, opOut.IsMontgomery, opOut.IsBatched = op0.IsNTT, op0.IsMontgomery, op0.IsBatched
	return nil
}

func (w *ScalarWorkspace) MulThenAddRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows, levelsConsumed int, scalar *bignum.Complex) error {
	level, err := w.validateUnary(ringQ, op0, opOut, rows)
	if err != nil {
		return err
	}
	if scalar == nil {
		return errors.New("Fast scalar MulThenAdd scalar cannot be nil")
	}
	if op0.Degree() != opOut.Degree() {
		return errors.New("Fast scalar MulThenAdd requires matching input/output degrees")
	}
	if len(op0.Value) > len(w.scratch) {
		return fmt.Errorf("Fast scalar MulThenAdd degree %d exceeds shared scratch capacity %d", op0.Degree(), len(w.scratch)-1)
	}

	scale := rlwe.Scale{}
	cmp := op0.Scale.Cmp(opOut.Scale)
	var promote *big.Int
	if cmp == 0 {
		if scalar.IsInt() {
			scale = rlwe.NewScale(1)
		} else {
			scale, err = CoefficientScale(ringQ, level, levelsConsumed)
			if err != nil {
				return err
			}
			promote = scale.BigInt()
		}
	} else if cmp < 0 {
		scale = opOut.Scale.Div(op0.Scale)
	} else {
		promote = op0.Scale.Div(opOut.Scale).BigInt()
		if promote.Sign() <= 0 {
			return errors.New("invalid Fast MulThenAdd scale promotion ratio")
		}
		if scalar.IsInt() {
			scale = rlwe.NewScale(1)
		} else {
			scale, err = CoefficientScale(ringQ, level, levelsConsumed)
			if err != nil {
				return err
			}
			promote.Mul(promote, scale.BigInt())
		}
	}
	values, err := ScalarNTTRows(ringQ, scalar, &scale.Value, false, level, rows)
	if err != nil {
		return err
	}

	source := op0.Value
	if promote != nil && op0 == opOut {
		for component := range op0.Value {
			CopyPrefixRowsUnchecked(rows, op0.Value[component], w.scratch[component])
		}
		source = w.scratch[:len(op0.Value)]
	}
	if promote != nil {
		if err := MulIntegerRows(ringQ, opOut, opOut, rows, promote, true); err != nil {
			return err
		}
		if cmp > 0 {
			opOut.Scale = op0.Scale.Mul(scale)
		} else {
			opOut.Scale = opOut.Scale.Mul(scale)
		}
	}

	return w.MulThenAddEncodedRows(ringQ, source, opOut, rows, values)
}

func (w *ScalarWorkspace) MulThenAddEncodedRows(ringQ *ring.Ring, source []ring.Poly, opOut *rlwe.Ciphertext, rows int, values [MaxQPrefixWidth][2]uint64) error {
	if w == nil || w.n < 1 || ringQ == nil || opOut == nil || len(source) == 0 || opOut.MetaData == nil {
		return errors.New("Fast scalar MulThenAdd encoded operands cannot be nil")
	}
	if ringQ.N() != w.n || opOut.N() != w.n || !opOut.IsNTT || rows > w.maxRows {
		return errors.New("Fast scalar MulThenAdd encoded dimensions or representation are invalid")
	}
	level := opOut.Level()
	if err := ValidatePrefixBacking(ringQ, rows, source...); err != nil {
		return fmt.Errorf("Fast scalar MulThenAdd source: %w", err)
	}
	if err := ValidatePrefixRows(ringQ, level, rows, opOut.Value...); err != nil {
		return fmt.Errorf("Fast scalar MulThenAdd output: %w", err)
	}
	for component := range source {
		for limb := 0; limb < rows; limb++ {
			subring := ringQ.SubRings[limb]
			half := ringQ.N() >> 1
			scalar0 := ring.MForm(values[limb][0], subring.Modulus, subring.BRedConstant)
			scalar1 := ring.MForm(values[limb][1], subring.Modulus, subring.BRedConstant)
			subring.MulScalarMontgomeryThenAdd(source[component].Coeffs[limb][:half], scalar0, opOut.Value[component].Coeffs[limb][:half])
			subring.MulScalarMontgomeryThenAdd(source[component].Coeffs[limb][half:], scalar1, opOut.Value[component].Coeffs[limb][half:])
		}
	}
	return nil
}

func CoefficientScale(ringQ *ring.Ring, level, levelsConsumed int) (rlwe.Scale, error) {
	if ringQ == nil {
		return rlwe.Scale{}, errors.New("Fast scalar coefficient scale ring cannot be nil")
	}
	if levelsConsumed < 1 || level-levelsConsumed+1 < 0 || level > ringQ.Level() {
		return rlwe.Scale{}, fmt.Errorf("insufficient level %d for scalar encoding across %d rescaling levels", level, levelsConsumed)
	}
	scale := rlwe.NewScale(1)
	for i := 0; i < levelsConsumed; i++ {
		scale = scale.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
	}
	return scale, nil
}

func (w *ScalarWorkspace) validateUnary(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows int) (int, error) {
	if w == nil || w.n < 1 || w.maxRows < 1 {
		return 0, errors.New("Fast scalar workspace is not initialized")
	}
	if ringQ == nil || op0 == nil || opOut == nil {
		return 0, errors.New("Fast scalar ring and operands cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil || len(op0.Value) == 0 || len(opOut.Value) == 0 {
		return 0, errors.New("Fast scalar operands require metadata and polynomial storage")
	}
	if ringQ.N() != w.n || op0.N() != w.n || opOut.N() != w.n {
		return 0, errors.New("Fast scalar operand dimensions do not match workspace")
	}
	if !op0.IsNTT || !opOut.IsNTT || op0.IsNTT != opOut.IsNTT || op0.IsMontgomery != opOut.IsMontgomery {
		return 0, errors.New("Fast scalar operands require matching NTT/Montgomery representations")
	}
	level := min(op0.Level(), opOut.Level())
	if level < 0 {
		return 0, errors.New("Fast scalar arithmetic level cannot be negative")
	}
	if rows > w.maxRows {
		return 0, fmt.Errorf("Fast scalar row count %d exceeds workspace capacity %d", rows, w.maxRows)
	}
	if err := ValidatePrefixRows(ringQ, level, rows, op0.Value...); err != nil {
		return 0, fmt.Errorf("Fast scalar input: %w", err)
	}
	for component := range opOut.Value {
		if err := ValidatePrefixRows(ringQ, level, rows, opOut.Value[component]); err != nil {
			return 0, fmt.Errorf("Fast scalar output component %d: %w", component, err)
		}
	}
	return level, nil
}

func (w *ScalarWorkspace) resizeOutput(op0, opOut *rlwe.Ciphertext, level int) {
	ResizeCompactCiphertext(opOut, op0.Degree(), level, w.n)
}

// MulIntegerRows multiplies only the explicitly authorized Q-prefix rows by
// scalar, preserving the existing Level and Scale. When preserveHigherRows is
// true, rows above the prefix are neither read nor written; otherwise the
// output is resized to compact Q-prefix storage.
func MulIntegerRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, rows int, scalar *big.Int, preserveHigherRows bool) error {
	if ringQ == nil || op0 == nil || opOut == nil || scalar == nil {
		return errors.New("Fast integer multiplication ring, operands and scalar cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return errors.New("Fast integer multiplication metadata cannot be nil")
	}
	if op0.N() != ringQ.N() || opOut.N() != ringQ.N() {
		return errors.New("Fast integer multiplication dimensions do not match ring")
	}
	if !op0.IsNTT || !opOut.IsNTT || op0.IsMontgomery != opOut.IsMontgomery {
		return errors.New("Fast integer multiplication requires matching NTT/Montgomery representations")
	}
	level := min(op0.Level(), opOut.Level())
	if err := ValidatePrefixRows(ringQ, level, rows, op0.Value...); err != nil {
		return fmt.Errorf("input: %w", err)
	}
	if !preserveHigherRows {
		ResizeCompactCiphertext(opOut, op0.Degree(), level, ringQ.N())
	}
	if err := ValidatePrefixRows(ringQ, level, rows, opOut.Value...); err != nil {
		return fmt.Errorf("output: %w", err)
	}
	var scalarMontgomery [MaxQPrefixWidth]uint64
	for limb := 0; limb < rows; limb++ {
		subring := ringQ.SubRings[limb]
		scalarMod := new(big.Int).Mod(scalar, new(big.Int).SetUint64(subring.Modulus)).Uint64()
		scalarMontgomery[limb] = ring.MForm(scalarMod, subring.Modulus, subring.BRedConstant)
	}
	for component := range op0.Value {
		for limb := 0; limb < rows; limb++ {
			subring := ringQ.SubRings[limb]
			subring.MulScalarMontgomery(op0.Value[component].Coeffs[limb], scalarMontgomery[limb], opOut.Value[component].Coeffs[limb])
		}
	}
	*opOut.MetaData = *op0.MetaData
	return nil
}
