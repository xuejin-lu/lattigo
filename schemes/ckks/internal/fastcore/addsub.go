package fastcore

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// MaxQPrefixWidth is the current Q-prefix engineering cap shared by the
// explicit Fast evaluator and the ordinary public Fast adapter.
const MaxQPrefixWidth = 4

// AddSubCore is the import-neutral interface used by the public CKKS and
// explicit Fast evaluator facades for the same bounded ciphertext kernel.
type AddSubCore interface {
	ApplyRows(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, rows int, sub bool) error
}

// AddSubWorkspace owns no mutable scratch; it exists as a narrow injectable
// dispatch boundary so both evaluator facades can share and test this kernel.
type AddSubWorkspace struct{}

func NewAddSubWorkspace() *AddSubWorkspace { return &AddSubWorkspace{} }

func (w *AddSubWorkspace) ApplyRows(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, rows int, sub bool) error {
	if w == nil {
		return errors.New("Fast Add/Sub core is not initialized")
	}
	return AddSubCiphertextsRows(ringQ, op0, op1, opOut, rows, sub)
}

// QPrefixWidth returns the maintained actual-logical-Q row count at level.
func QPrefixWidth(level int) (int, error) {
	if level < 0 {
		return 0, fmt.Errorf("Q-prefix level must be non-negative: %d", level)
	}
	if level >= MaxQPrefixWidth-1 {
		return MaxQPrefixWidth, nil
	}
	return level + 1, nil
}

// AddSubPolynomialPair describes one component-wise modular add/subtract.
type AddSubPolynomialPair struct {
	Op0    ring.Poly
	Op1    ring.Poly
	Output ring.Poly
}

// AddSubPolynomialRows applies modular addition/subtraction to exactly rows
// for every pair. It validates every pair before writing any output.
func AddSubPolynomialRows(ringQ *ring.Ring, level, rows int, sub bool, pairs ...AddSubPolynomialPair) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if level < 0 || level > ringQ.Level() {
		return fmt.Errorf("Q-prefix level %d is outside ring level [0,%d]", level, ringQ.Level())
	}
	width, err := QPrefixWidth(level)
	if err != nil {
		return err
	}
	if rows < 1 || rows > level+1 || rows > width || rows > len(ringQ.SubRings) {
		return fmt.Errorf("Q-prefix row count %d must be in [1,%d] at level %d", rows, min(level+1, width), level)
	}
	if len(pairs) == 0 {
		return errors.New("Fast Add/Sub requires at least one polynomial pair")
	}
	for row := 0; row < rows; row++ {
		if ringQ.SubRings[row] == nil {
			return fmt.Errorf("q%d subring is unavailable", row)
		}
	}
	for i, pair := range pairs {
		for _, poly := range []ring.Poly{pair.Op0, pair.Op1, pair.Output} {
			if len(poly.Coeffs) < level+1 {
				return fmt.Errorf("polynomial pair %d has %d logical rows, requires level %d", i, len(poly.Coeffs), level)
			}
			for row := 0; row < rows; row++ {
				if len(poly.Coeffs[row]) != ringQ.N() {
					return fmt.Errorf("polynomial pair %d q%d backing length %d does not match N=%d", i, row, len(poly.Coeffs[row]), ringQ.N())
				}
			}
		}
	}
	for _, pair := range pairs {
		for row := 0; row < rows; row++ {
			if sub {
				ringQ.SubRings[row].Sub(pair.Op0.Coeffs[row], pair.Op1.Coeffs[row], pair.Output.Coeffs[row])
			} else {
				ringQ.SubRings[row].Add(pair.Op0.Coeffs[row], pair.Op1.Coeffs[row], pair.Output.Coeffs[row])
			}
		}
	}
	return nil
}

// AddSubCiphertextsRows applies the existing Fast ciphertext add/subtract
// contract to exactly rows Q-prefix limbs and compacts the output backing.
func AddSubCiphertextsRows(ringQ *ring.Ring, op0, op1, opOut *rlwe.Ciphertext, rows int, sub bool) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if op0 == nil || op1 == nil || opOut == nil {
		return errors.New("op0, op1 and opOut cannot be nil")
	}
	if len(op0.Value) == 0 || len(op1.Value) == 0 || len(opOut.Value) == 0 {
		return errors.New("ciphertexts must contain polynomial storage")
	}
	if op0.MetaData == nil || op1.MetaData == nil || opOut.MetaData == nil {
		return errors.New("op0, op1 and opOut metadata cannot be nil")
	}
	if ringQ.Level() < 1 || op0.Level() < 1 || op1.Level() < 1 || opOut.Level() < 1 {
		return errors.New("FastAdd/FastSub requires q0 and q1")
	}
	if !op0.Scale.Equal(op1.Scale) {
		return errors.New("FastAdd/FastSub requires equal operand scales")
	}
	if op0.IsNTT != op1.IsNTT {
		return errors.New("FastAdd/FastSub requires equal IsNTT domains")
	}
	if op0.IsMontgomery != op1.IsMontgomery {
		return errors.New("FastAdd/FastSub requires equal Montgomery representations")
	}

	level := utils.Min(op0.Level(), op1.Level())
	level = utils.Min(level, opOut.Level())
	if level < 1 {
		return errors.New("FastAdd/FastSub output level must contain q0 and q1")
	}
	maxDegree := utils.Max(op0.Degree(), op1.Degree())
	minDegree := utils.Min(op0.Degree(), op1.Degree())
	validate := func(name string, ct *rlwe.Ciphertext, degree int) error {
		if ct.N() != opOut.N() || ct.N() != ringQ.N() {
			return fmt.Errorf("%s dimension does not match output", name)
		}
		for component := 0; component <= degree; component++ {
			if component >= len(ct.Value) {
				return fmt.Errorf("%s is missing component %d", name, component)
			}
			if err := ValidatePrefixRows(ringQ, level, rows, ct.Value[component]); err != nil {
				return fmt.Errorf("%s component %d: %w", name, component, err)
			}
		}
		return nil
	}
	if err := validate("op0", op0, op0.Degree()); err != nil {
		return err
	}
	if err := validate("op1", op1, op1.Degree()); err != nil {
		return err
	}
	if opOut.N() == 0 || opOut.N() != ringQ.N() {
		return errors.New("output ciphertext has invalid storage")
	}
	if err := validate("opOut", opOut, min(opOut.Degree(), maxDegree)); err != nil {
		return err
	}

	ResizeCompactCiphertext(opOut, maxDegree, level, ringQ.N())
	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale
	opOut.IsNTT = op0.IsNTT
	opOut.IsMontgomery = op0.IsMontgomery
	opOut.IsBatched = op0.IsBatched
	opOut.LogDimensions.Rows = utils.Max(op0.LogDimensions.Rows, op1.LogDimensions.Rows)
	opOut.LogDimensions.Cols = utils.Max(op0.LogDimensions.Cols, op1.LogDimensions.Cols)

	pairs := make([]AddSubPolynomialPair, minDegree+1)
	for component := 0; component <= minDegree; component++ {
		pairs[component] = AddSubPolynomialPair{Op0: op0.Value[component], Op1: op1.Value[component], Output: opOut.Value[component]}
	}
	if err := AddSubPolynomialRows(ringQ, level, rows, sub, pairs...); err != nil {
		return err
	}

	if op0.Degree() > minDegree && opOut != op0 {
		for component := minDegree + 1; component <= op0.Degree(); component++ {
			CopyPrefixRowsUnchecked(rows, op0.Value[component], opOut.Value[component])
		}
	} else if op1.Degree() > minDegree {
		for component := minDegree + 1; component <= op1.Degree(); component++ {
			if sub {
				if err := NegatePrefixRows(ringQ, level, rows, op1.Value[component], opOut.Value[component]); err != nil {
					return err
				}
			} else if opOut != op1 {
				CopyPrefixRowsUnchecked(rows, op1.Value[component], opOut.Value[component])
			}
		}
	}
	return nil
}

// NewCompactCiphertext allocates logical-level ciphertext metadata with only
// the current capped Q-prefix rows backed by N coefficients.
func NewCompactCiphertext(params rlwe.ParameterProvider, degree, level int) *rlwe.Ciphertext {
	n := params.GetRLWEParameters().N()
	prefixWidth, err := QPrefixWidth(level)
	if err != nil {
		panic(err)
	}
	polys := make([]ring.Poly, degree+1)
	for i := range polys {
		polys[i] = compactPoly(n, level, prefixWidth)
	}
	ct, err := rlwe.NewCiphertextAtLevelFromPoly(level, polys)
	if err != nil {
		panic(err)
	}
	ct.IsNTT = params.GetRLWEParameters().NTTFlag()
	return ct
}

// ResizeCompactCiphertext changes logical level/degree without materializing
// dormant rows above the current Q-prefix.
func ResizeCompactCiphertext(ct *rlwe.Ciphertext, degree, level, n int) {
	if ct == nil || level < 0 || degree < 0 {
		panic("invalid Fast ciphertext resize")
	}
	prefixWidth, err := QPrefixWidth(level)
	if err != nil {
		panic(err)
	}
	for i := range ct.Value {
		coeffs := ct.Value[i].Coeffs
		if level < len(coeffs)-1 {
			ct.Value[i].Coeffs = coeffs[:level+1]
		} else if level > len(coeffs)-1 {
			grown := make([][]uint64, level+1)
			copy(grown, coeffs)
			ct.Value[i].Coeffs = grown
		}
		ensureQPrefixRows(&ct.Value[i], level, n, prefixWidth)
	}
	if degree < len(ct.Value)-1 {
		ct.Value = ct.Value[:degree+1]
	} else {
		for len(ct.Value) < degree+1 {
			ct.Value = append(ct.Value, compactPoly(n, level, prefixWidth))
		}
	}
}

func compactPoly(n, level, prefixWidth int) ring.Poly {
	coeffs := make([][]uint64, level+1)
	for row := 0; row < len(coeffs) && row < prefixWidth; row++ {
		coeffs[row] = make([]uint64, n)
	}
	return ring.Poly{Coeffs: coeffs}
}

func ensureQPrefixRows(poly *ring.Poly, level, n, prefixWidth int) {
	for row := range poly.Coeffs {
		if row < prefixWidth {
			if len(poly.Coeffs[row]) != n {
				poly.Coeffs[row] = make([]uint64, n)
			}
		} else {
			poly.Coeffs[row] = nil
		}
	}
}

// ValidatePrefixRows checks logical-level and physical backing for an
// explicitly bounded Q-prefix operation.
func ValidatePrefixRows(ringQ *ring.Ring, level, rows int, polys ...ring.Poly) error {
	if ringQ == nil {
		return errors.New("Q-prefix ring cannot be nil")
	}
	if level < 0 || level > ringQ.Level() {
		return fmt.Errorf("Q-prefix level %d is outside ring level [0,%d]", level, ringQ.Level())
	}
	width, err := QPrefixWidth(level)
	if err != nil {
		return err
	}
	if rows < 1 || rows > level+1 || rows > width {
		return fmt.Errorf("Q-prefix row count %d must be in [1,%d] at level %d", rows, min(level+1, width), level)
	}
	if rows > len(ringQ.SubRings) {
		return fmt.Errorf("Q-prefix row count %d exceeds available q subrings %d", rows, len(ringQ.SubRings))
	}
	for row := 0; row < rows; row++ {
		if ringQ.SubRings[row] == nil {
			return fmt.Errorf("q%d subring is unavailable", row)
		}
	}
	for i, poly := range polys {
		if len(poly.Coeffs) == 0 || len(poly.Coeffs) < level+1 {
			return fmt.Errorf("polynomial %d logical level is below requested level %d", i, level)
		}
	}
	return ValidatePrefixBacking(ringQ, rows, polys...)
}

// ValidatePrefixBacking checks physical row availability without inferring
// authority from allocated rows above the explicit request.
func ValidatePrefixBacking(ringQ *ring.Ring, rows int, polys ...ring.Poly) error {
	if ringQ == nil {
		return errors.New("Q-prefix ring cannot be nil")
	}
	if rows < 1 || rows > MaxQPrefixWidth || rows > len(ringQ.SubRings) {
		return fmt.Errorf("Q-prefix backing row count %d is unavailable", rows)
	}
	for row := 0; row < rows; row++ {
		if ringQ.SubRings[row] == nil {
			return fmt.Errorf("q%d subring is unavailable", row)
		}
	}
	for i, poly := range polys {
		if len(poly.Coeffs) == 0 {
			return fmt.Errorf("polynomial %d has no coefficient rows", i)
		}
		if poly.N() != ringQ.N() {
			return fmt.Errorf("polynomial %d dimension %d does not match ring dimension %d", i, poly.N(), ringQ.N())
		}
		if len(poly.Coeffs) < rows {
			return fmt.Errorf("polynomial %d has %d rows, requested %d", i, len(poly.Coeffs), rows)
		}
		for row := 0; row < rows; row++ {
			if len(poly.Coeffs[row]) != ringQ.N() {
				return fmt.Errorf("polynomial %d q%d backing length %d does not match N=%d", i, row, len(poly.Coeffs[row]), ringQ.N())
			}
		}
	}
	return nil
}

// CopyPrefixRowsUnchecked copies exactly rows already-validated polynomial
// limbs; rows above the request are not read or written.
func CopyPrefixRowsUnchecked(rows int, src, dst ring.Poly) {
	for row := 0; row < rows; row++ {
		copy(dst.Coeffs[row], src.Coeffs[row])
	}
}

// CopyPrefixRows validates and copies exactly rows; rows above the request
// are neither read nor written.
func CopyPrefixRows(ringQ *ring.Ring, level, rows int, src, dst ring.Poly) error {
	if err := ValidatePrefixRows(ringQ, level, rows, src, dst); err != nil {
		return err
	}
	CopyPrefixRowsUnchecked(rows, src, dst)
	return nil
}

// NegatePrefixRows negates exactly rows after validating all backing.
func NegatePrefixRows(ringQ *ring.Ring, level, rows int, src, dst ring.Poly) error {
	if err := ValidatePrefixRows(ringQ, level, rows, src, dst); err != nil {
		return err
	}
	for row := 0; row < rows; row++ {
		ringQ.SubRings[row].Neg(src.Coeffs[row], dst.Coeffs[row])
	}
	return nil
}
