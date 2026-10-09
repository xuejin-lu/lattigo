package fastcore

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// MulCore is the import-neutral boundary shared by the ordinary public CKKS
// evaluator's zero-secret path and the explicit Fast evaluator.
type MulCore interface {
	ApplyRows(ringQ *ring.Ring, op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], opOut *rlwe.Ciphertext, rows int, relin bool) error
}

// MulWorkspace owns the three temporary polynomials used by the existing
// pointwise NTT product kernel. A workspace is single-stream, like its
// evaluator owner.
type MulWorkspace struct {
	n       int
	maxRows int
	scratch [3]ring.Poly
}

func NewMulWorkspace(ringQ *ring.Ring, maxRows int) *MulWorkspace {
	w := &MulWorkspace{maxRows: maxRows}
	if ringQ == nil || maxRows < 1 || maxRows > MaxQPrefixWidth {
		return w
	}
	w.n = ringQ.N()
	for i := range w.scratch {
		w.scratch[i] = ring.NewPoly(w.n, maxRows-1)
	}
	return w
}

func (w *MulWorkspace) ApplyRows(ringQ *ring.Ring, op0 *rlwe.Ciphertext, op1 *rlwe.Element[ring.Poly], opOut *rlwe.Ciphertext, rows int, relin bool) error {
	if w == nil || w.n < 1 || w.maxRows < 1 {
		return errors.New("Fast Mul workspace is not initialized")
	}
	if ringQ == nil || op0 == nil || op1 == nil || opOut == nil {
		return errors.New("Fast Mul ring and operands cannot be nil")
	}
	if op0.MetaData == nil || op1.MetaData == nil || opOut.MetaData == nil {
		return errors.New("Fast Mul operand metadata cannot be nil")
	}
	if len(op0.Value) == 0 || len(op1.Value) == 0 {
		return errors.New("Fast Mul inputs must contain polynomial storage")
	}
	if op0.N() != w.n || op1.N() != w.n || opOut.N() != w.n || ringQ.N() != w.n {
		return errors.New("Fast Mul operand dimensions do not match workspace")
	}
	if op0.Degree() > 1 || op1.Degree() > 1 {
		return errors.New("Fast Mul requires operands of degree at most one")
	}
	if !op0.IsNTT || !op1.IsNTT || op0.IsNTT != op1.IsNTT || op0.IsMontgomery != op1.IsMontgomery {
		return errors.New("Fast Mul operands require matching NTT/Montgomery representations")
	}
	if op0.IsBatched != op1.IsBatched {
		return errors.New("Fast Mul operands require matching batching metadata")
	}

	level := utils.Min(op0.Level(), utils.Min(op1.Level(), opOut.Level()))
	if err := ValidatePrefixRows(ringQ, level, rows, op0.Value...); err != nil {
		return fmt.Errorf("Fast Mul op0: %w", err)
	}
	if err := ValidatePrefixRows(ringQ, level, rows, op1.Value...); err != nil {
		return fmt.Errorf("Fast Mul op1: %w", err)
	}
	if rows > w.maxRows || rows > len(w.scratch[0].Coeffs) {
		return fmt.Errorf("Fast Mul row count %d exceeds workspace capacity %d", rows, w.maxRows)
	}
	for i := range w.scratch {
		if err := ValidatePrefixBacking(ringQ, rows, w.scratch[i]); err != nil {
			return fmt.Errorf("Fast Mul scratch %d: %w", i, err)
		}
	}
	for component, polys := range [][]ring.Poly{op0.Value, op1.Value} {
		for d, poly := range polys {
			for row := 0; row < rows; row++ {
				modulus := ringQ.SubRings[row].Modulus
				for coefficient, value := range poly.Coeffs[row] {
					if value >= modulus {
						return fmt.Errorf("Fast Mul input set %d component %d q%d coefficient %d is not a canonical residue", component, d, row, coefficient)
					}
				}
			}
		}
	}
	if len(opOut.Value) == 0 {
		return errors.New("Fast Mul output has no polynomial storage")
	}
	for d := range opOut.Value {
		if err := ValidatePrefixRows(ringQ, level, rows, opOut.Value[d]); err != nil {
			return fmt.Errorf("Fast Mul output component %d: %w", d, err)
		}
	}

	degree := op0.Degree() + op1.Degree()
	if relin && degree == 2 {
		degree = 1
	}
	if op0.Degree() == 1 && op1.Degree() == 1 {
		t0, t1, t2 := w.scratch[0], w.scratch[1], w.scratch[2]
		if err := mulPrefixRows(ringQ, rows, op0.Value[0], op1.Value[0], t0, op0.IsMontgomery, false); err != nil {
			return err
		}
		if op0.El() == op1 {
			if err := mulPrefixRows(ringQ, rows, op0.Value[0], op0.Value[1], t1, op0.IsMontgomery, false); err != nil {
				return err
			}
			for row := 0; row < rows; row++ {
				ringQ.SubRings[row].Add(t1.Coeffs[row], t1.Coeffs[row], t1.Coeffs[row])
			}
		} else {
			if err := mulPrefixRows(ringQ, rows, op0.Value[0], op1.Value[1], t1, op0.IsMontgomery, false); err != nil {
				return err
			}
			if err := mulPrefixRows(ringQ, rows, op0.Value[1], op1.Value[0], t1, op0.IsMontgomery, true); err != nil {
				return err
			}
		}
		if !relin {
			if err := mulPrefixRows(ringQ, rows, op0.Value[1], op1.Value[1], t2, op0.IsMontgomery, false); err != nil {
				return err
			}
		}
		ResizeCompactCiphertext(opOut, degree, level, w.n)
		CopyPrefixRowsUnchecked(rows, t0, opOut.Value[0])
		CopyPrefixRowsUnchecked(rows, t1, opOut.Value[1])
		if !relin {
			CopyPrefixRowsUnchecked(rows, t2, opOut.Value[2])
		}
	} else {
		var ct []ring.Poly
		var pt ring.Poly
		if op0.Degree() == 0 {
			pt, ct = op0.Value[0], op1.Value
		} else {
			pt, ct = op1.Value[0], op0.Value
		}
		for d := range ct {
			if err := mulPrefixRows(ringQ, rows, pt, ct[d], w.scratch[d], op0.IsMontgomery, false); err != nil {
				return err
			}
		}
		ResizeCompactCiphertext(opOut, degree, level, w.n)
		for d := range ct {
			CopyPrefixRowsUnchecked(rows, w.scratch[d], opOut.Value[d])
		}
	}

	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale.Mul(op1.Scale)
	opOut.IsNTT = op0.IsNTT
	opOut.IsMontgomery = op0.IsMontgomery
	opOut.IsBatched = op0.IsBatched
	opOut.LogDimensions.Rows = utils.Max(op0.LogDimensions.Rows, op1.LogDimensions.Rows)
	opOut.LogDimensions.Cols = utils.Max(op0.LogDimensions.Cols, op1.LogDimensions.Cols)
	return nil
}

func mulPrefixRows(ringQ *ring.Ring, rows int, a, b, out ring.Poly, montgomery, add bool) error {
	if err := ValidatePrefixBacking(ringQ, rows, a, b, out); err != nil {
		return err
	}
	for row := 0; row < rows; row++ {
		subring := ringQ.SubRings[row]
		switch {
		case montgomery && add:
			subring.MulCoeffsMontgomeryThenAdd(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		case montgomery:
			subring.MulCoeffsMontgomery(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		case add:
			subring.MulCoeffsBarrettThenAdd(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		default:
			subring.MulCoeffsBarrett(a.Coeffs[row], b.Coeffs[row], out.Coeffs[row])
		}
	}
	return nil
}
