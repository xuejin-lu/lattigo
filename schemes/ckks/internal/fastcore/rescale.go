package fastcore

import (
	"errors"
	"fmt"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/internal/fastdiag"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// Uint192 is the fixed-width representation used by the existing bounded
// Q-prefix centered-CRT rescale kernel.
type Uint192 struct{ Lo, Mid, Hi uint64 }

// QPrefixCapacityError reports the first output component that does not fit
// the strict centered interval at a logical level transition.
type QPrefixCapacityError struct {
	Level         int
	Component     int
	Bound         *big.Int
	PrefixProduct *big.Int
}

func (err *QPrefixCapacityError) Error() string {
	if err == nil {
		return "Q-prefix capacity exceeded"
	}
	bound := "<nil>"
	if err.Bound != nil {
		bound = new(big.Int).Lsh(new(big.Int).Set(err.Bound), 1).String()
	}
	product := "<nil>"
	if err.PrefixProduct != nil {
		product = err.PrefixProduct.String()
	}
	return fmt.Sprintf("Q-prefix capacity exceeded at Level %d component %d: strict 2B < S_Q failed (2B=%s, S_Q=%s)", err.Level, err.Component, bound, product)
}

// RescaleCore exposes the existing bounded centered-CRT rescale implementation
// to the public CKKS and explicit Fast evaluator facades.
type RescaleCore interface {
	ApplyRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, nbRescales, sourceRows int) error
	ApplyToRows(ringQ *ring.Ring, op0 *rlwe.Ciphertext, minScale rlwe.Scale, sourceRows int, opOut *rlwe.Ciphertext) error
}

// RescaleWorkspace holds evaluator-local scratch for transactional Q-prefix
// Rescale. A workspace is single-stream, matching the owning evaluator.
type RescaleWorkspace struct {
	n        int
	maxRows  int
	coeff    ring.Poly
	staged   []ring.Poly
	diagPar  uint64
	q        [MaxQPrefixWidth]uint64
	inverse  [MaxQPrefixWidth]uint64
	modulus  [MaxQPrefixWidth]Uint192
	half     [MaxQPrefixWidth]Uint192
	initErr  error
	widthErr [MaxQPrefixWidth]error
}

func NewRescaleWorkspace(ringQ *ring.Ring, maxRows int) *RescaleWorkspace {
	w := &RescaleWorkspace{maxRows: maxRows}
	if ringQ == nil || maxRows < 1 || maxRows > MaxQPrefixWidth || maxRows > ringQ.Level()+1 {
		w.initErr = fmt.Errorf("invalid Q-prefix Rescale workspace width %d", maxRows)
		return w
	}
	w.n = ringQ.N()
	w.coeff = ring.NewPoly(w.n, maxRows-1)
	w.staged = make([]ring.Poly, 2)
	for i := range w.staged {
		w.staged[i] = ring.NewPoly(w.n, maxRows-1)
	}
	product := Uint192{Lo: 1}
	productValid := true
	for row := 0; row < maxRows; row++ {
		if ringQ.SubRings[row] == nil {
			w.widthErr[row] = fmt.Errorf("Fast Rescale q%d subring is unavailable", row)
			productValid = false
			continue
		}
		q := ringQ.SubRings[row].Modulus
		if q < 3 || q&1 == 0 || q > uint64(^uint64(0)>>1) {
			w.widthErr[row] = fmt.Errorf("unsupported Fast Rescale modulus q%d=%d", row, q)
			productValid = false
			continue
		}
		w.q[row] = q
		if !productValid {
			w.widthErr[row] = fmt.Errorf("Fast Rescale Q-prefix through q%d is unavailable", row)
			continue
		}
		var overflow bool
		product, overflow = mul192By64(product, q)
		if overflow {
			w.widthErr[row] = fmt.Errorf("Fast Rescale Q-prefix product through q%d exceeds 192 bits", row)
			productValid = false
			continue
		}
		w.modulus[row] = product
		w.half[row] = half192(product)
		if row > 0 {
			inverse, ok := inverseMod(mod192By64(w.modulus[row-1], q), q)
			if !ok {
				w.widthErr[row] = fmt.Errorf("Fast Rescale Q-prefix through q%d is not pairwise coprime", row)
				productValid = false
				continue
			}
			w.inverse[row] = inverse
		}
	}
	return w
}

func (w *RescaleWorkspace) validateWidth(rows int) error {
	if w == nil {
		return errors.New("Fast Rescale workspace cannot be nil")
	}
	if w.initErr != nil {
		return w.initErr
	}
	if rows < 1 || rows > w.maxRows || rows > MaxQPrefixWidth {
		return fmt.Errorf("invalid Fast Rescale Q-prefix width %d", rows)
	}
	return w.widthErr[rows-1]
}

func (w *RescaleWorkspace) ensureStaged(ringQ *ring.Ring, components int) error {
	if w == nil || ringQ == nil {
		return errors.New("Fast Rescale workspace and ring cannot be nil")
	}
	if components < 1 {
		return errors.New("Fast Rescale ciphertext has no components")
	}
	for len(w.staged) < components {
		w.staged = append(w.staged, ring.NewPoly(ringQ.N(), w.maxRows-1))
	}
	for component := range w.staged {
		if len(w.staged[component].Coeffs) < w.maxRows || w.staged[component].N() != ringQ.N() {
			w.staged[component] = ring.NewPoly(ringQ.N(), w.maxRows-1)
		}
	}
	return nil
}

// ApplyRows applies Standard CKKS Rescale semantics using the logical top-q
// divisor and exactly sourceRows authoritative source residues. It stages all
// components and checks every intermediate centered-capacity boundary before
// committing any caller-visible output.
func (w *RescaleWorkspace) ApplyRows(ringQ *ring.Ring, op0, opOut *rlwe.Ciphertext, nbRescales, sourceRows int) error {
	var wholeSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		wholeSpan = fastdiag.Begin(fastdiag.Rescale, "rescale", 0, fastdiag.Input(op0, sourceRows).Merge(fastdiag.InPlace(op0 == opOut)).Merge(fastdiag.RepetitionCount(nbRescales)))
	}
	if w == nil || ringQ == nil || op0 == nil || opOut == nil {
		return errors.New("Fast Rescale workspace, ring and ciphertexts cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return errors.New("Fast Rescale metadata cannot be nil")
	}
	if nbRescales <= 0 || op0.Level() < nbRescales {
		return errors.New("invalid number of Fast rescale levels")
	}
	if op0.N() != w.n || opOut.N() != w.n || ringQ.N() != w.n {
		return errors.New("ciphertext dimensions do not match Fast Rescale workspace")
	}
	if !op0.IsNTT {
		return errors.New("Fast Rescale requires NTT-domain ciphertexts")
	}
	if len(op0.Value) == 0 || op0.Level() < 1 {
		return errors.New("Fast Rescale requires a non-zero logical Level and at least one component")
	}
	if ringQ.Level() < op0.Level() {
		return fmt.Errorf("Fast Rescale input Level %d exceeds configured ring Level %d", op0.Level(), ringQ.Level())
	}
	maxSourceRows, err := QPrefixWidth(op0.Level())
	if err != nil {
		return err
	}
	if sourceRows < 1 || sourceRows > maxSourceRows || sourceRows > w.maxRows {
		return fmt.Errorf("explicit Fast Rescale source width %d is outside [1,%d] at Level %d", sourceRows, min(maxSourceRows, w.maxRows), op0.Level())
	}
	if err := w.validateWidth(sourceRows); err != nil {
		return err
	}
	targetLevel := op0.Level() - nbRescales
	maxTargetRows, err := QPrefixWidth(targetLevel)
	if err != nil {
		return err
	}
	targetRows := min(sourceRows, maxTargetRows)
	if err := w.ensureStaged(ringQ, len(op0.Value)); err != nil {
		return err
	}
	if len(op0.Value) == 0 || sourceRows > len(w.coeff.Coeffs) || targetRows > len(w.staged[0].Coeffs) {
		return fmt.Errorf("Fast Rescale workspace does not cover source/target widths %d/%d", sourceRows, targetRows)
	}
	for component := range op0.Value {
		if err := ValidatePrefixRows(ringQ, op0.Level(), sourceRows, op0.Value[component]); err != nil {
			return fmt.Errorf("Fast Rescale component %d: %w", component, err)
		}
		for row := 0; row < sourceRows; row++ {
			modulus := ringQ.SubRings[row].Modulus
			for coefficient, value := range op0.Value[component].Coeffs[row] {
				if value >= modulus {
					return fmt.Errorf("Fast Rescale component %d q%d coefficient %d is not a canonical residue", component, row, coefficient)
				}
			}
		}
	}
	for step := 0; step < nbRescales; step++ {
		logicalLevel := op0.Level() - step
		if logicalLevel >= len(ringQ.SubRings) || ringQ.SubRings[logicalLevel] == nil || ringQ.SubRings[logicalLevel].Modulus < 2 {
			return fmt.Errorf("Fast Rescale logical divisor q%d is unavailable", logicalLevel)
		}
	}

	var preflightSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		preflightSpan = fastdiag.Begin(fastdiag.Rescale, "preflight", wholeSpan.Sequence(), fastdiag.Input(op0, sourceRows).Merge(fastdiag.InPlace(op0 == opOut)).Merge(fastdiag.RepetitionCount(len(op0.Value))))
		w.diagPar = preflightSpan.Sequence()
	}
	for component := range op0.Value {
		if err := w.rescaleComponent(ringQ, op0.Value[component], &w.staged[component], op0.Level(), nbRescales, sourceRows, targetRows, component, op0.IsMontgomery, op0 == opOut); err != nil {
			if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
				preflightSpan.End(fastdiag.Fields{})
			}
			return err
		}
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		preflightSpan.End(fastdiag.Output(op0, sourceRows))
	}

	var materializationSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		materializationSpan = fastdiag.Begin(fastdiag.Rescale, "materialization", wholeSpan.Sequence(), fastdiag.Input(op0, sourceRows).Merge(fastdiag.InPlace(op0 == opOut)).Merge(fastdiag.RepetitionCount(len(op0.Value))))
	}
	if opOut != op0 {
		ResizeCompactCiphertext(opOut, op0.Degree(), targetLevel, w.n)
	}
	for component := range op0.Value {
		var restoreSpan fastdiag.Span
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			restoreSpan = fastdiag.Begin(fastdiag.Rescale, "ntt_montgomery_restore", materializationSpan.Sequence(), rescaleDiagFields(op0.Level(), targetLevel, sourceRows, targetRows, fmt.Sprintf("c%d", component), op0 == opOut))
		}
		for row := 0; row < targetRows; row++ {
			ringQ.SubRings[row].NTT(w.staged[component].Coeffs[row], opOut.Value[component].Coeffs[row])
			if op0.IsMontgomery {
				ringQ.SubRings[row].MForm(opOut.Value[component].Coeffs[row], opOut.Value[component].Coeffs[row])
			}
		}
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			restoreSpan.End(fastdiag.Fields{})
		}
	}
	*opOut.MetaData = *op0.MetaData
	scale := op0.Scale
	for step := 0; step < nbRescales; step++ {
		scale = scale.Div(rlwe.NewScale(ringQ.SubRings[op0.Level()-step].Modulus))
	}
	opOut.Scale = scale
	opOut.IsNTT = op0.IsNTT
	opOut.IsMontgomery = op0.IsMontgomery
	opOut.IsBatched = op0.IsBatched
	if opOut == op0 {
		ResizeCompactCiphertext(opOut, op0.Degree(), targetLevel, w.n)
	}
	retainQPrefixRows(opOut, targetRows)
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		materializationSpan.End(fastdiag.Output(opOut, targetRows))
		wholeSpan.End(fastdiag.Output(opOut, targetRows))
	}
	return nil
}

func (w *RescaleWorkspace) ApplyToRows(ringQ *ring.Ring, op0 *rlwe.Ciphertext, minScale rlwe.Scale, sourceRows int, opOut *rlwe.Ciphertext) error {
	if w == nil || ringQ == nil || op0 == nil || opOut == nil {
		return errors.New("Fast RescaleTo workspace, ring and ciphertexts cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return errors.New("Fast RescaleTo metadata cannot be nil")
	}
	if minScale.Cmp(rlwe.NewScale(0)) != 1 || op0.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return errors.New("cannot RescaleTo: scales must be positive")
	}
	if op0.Level() < 1 {
		return errors.New("cannot RescaleTo: input Ciphertext already at level 0")
	}
	if op0.N() != w.n || opOut.N() != w.n || ringQ.N() != w.n || !op0.IsNTT {
		return errors.New("Fast RescaleTo dimensions or NTT domain are invalid")
	}
	maxSourceRows, err := QPrefixWidth(op0.Level())
	if err != nil {
		return err
	}
	if sourceRows < 1 || sourceRows > maxSourceRows || sourceRows > w.maxRows {
		return fmt.Errorf("explicit Fast RescaleTo source width %d is outside [1,%d] at Level %d", sourceRows, min(maxSourceRows, w.maxRows), op0.Level())
	}
	if err := w.validateWidth(sourceRows); err != nil {
		return err
	}
	if ringQ.Level() < op0.Level() {
		return fmt.Errorf("Fast RescaleTo input Level %d exceeds configured ring Level %d", op0.Level(), ringQ.Level())
	}
	for component := range op0.Value {
		if err := ValidatePrefixRows(ringQ, op0.Level(), sourceRows, op0.Value[component]); err != nil {
			return fmt.Errorf("Fast RescaleTo component %d: %w", component, err)
		}
		for row := 0; row < sourceRows; row++ {
			modulus := ringQ.SubRings[row].Modulus
			for coefficient, value := range op0.Value[component].Coeffs[row] {
				if value >= modulus {
					return fmt.Errorf("Fast RescaleTo component %d q%d coefficient %d is not a canonical residue", component, row, coefficient)
				}
			}
		}
	}
	threshold := minScale.Div(rlwe.NewScale(2))
	scale := op0.Scale
	newLevel, nbRescales := op0.Level(), 0
	for newLevel > 0 {
		if newLevel >= len(ringQ.SubRings) || ringQ.SubRings[newLevel] == nil || ringQ.SubRings[newLevel].Modulus < 2 {
			return fmt.Errorf("Fast RescaleTo logical divisor q%d is unavailable", newLevel)
		}
		candidate := scale.Div(rlwe.NewScale(ringQ.SubRings[newLevel].Modulus))
		if candidate.Cmp(threshold) == -1 {
			break
		}
		scale = candidate
		newLevel--
		nbRescales++
	}
	if nbRescales == 0 {
		if op0 != opOut {
			if err := w.ensureStaged(ringQ, len(op0.Value)); err != nil {
				return err
			}
			ResizeCompactCiphertext(opOut, op0.Degree(), op0.Level(), w.n)
			for component := range op0.Value {
				CopyPrefixRowsUnchecked(sourceRows, op0.Value[component], opOut.Value[component])
			}
			*opOut.MetaData = *op0.MetaData
			opOut.Scale = op0.Scale
			opOut.IsNTT, opOut.IsMontgomery, opOut.IsBatched = op0.IsNTT, op0.IsMontgomery, op0.IsBatched
		}
		// Explicit row-authority callers may intentionally provide fewer rows
		// than the ordinary Q-prefix policy. Keep those rows authoritative in
		// the output instead of exposing newly allocated zero rows as valid.
		retainQPrefixRows(opOut, sourceRows)
		return nil
	}
	return w.ApplyRows(ringQ, op0, opOut, nbRescales, sourceRows)
}

// retainQPrefixRows keeps the ciphertext's logical Level while marking every
// row outside the caller-authorized prefix as unavailable backing.
func retainQPrefixRows(ct *rlwe.Ciphertext, rows int) {
	for component := range ct.Value {
		for row := rows; row < len(ct.Value[component].Coeffs); row++ {
			ct.Value[component].Coeffs[row] = nil
		}
	}
}

func (w *RescaleWorkspace) rescaleComponent(ringQ *ring.Ring, src ring.Poly, staged *ring.Poly, sourceLevel, nbRescales, sourceRows, targetRows, component int, montgomery, inPlace bool) error {
	if staged == nil {
		return errors.New("Fast Rescale staging polynomial cannot be nil")
	}
	if err := ValidatePrefixBacking(ringQ, targetRows, *staged); err != nil {
		return fmt.Errorf("Fast Rescale component %d staging: %w", component, err)
	}
	fields := rescaleDiagFields(sourceLevel, sourceLevel-nbRescales, sourceRows, targetRows, fmt.Sprintf("c%d", component), inPlace)
	var phase fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		phase = fastdiag.Begin(fastdiag.Rescale, "prefix_to_coefficient", w.diagPar, fields)
	}
	if err := prefixToCoefficientRows(ringQ, src, sourceRows, true, montgomery, w.coeff); err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			phase.End(fastdiag.Fields{})
		}
		return err
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		phase.End(fastdiag.Fields{})
	}
	var reconstructNS, materializeNS int64
	var loop fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		loop = fastdiag.Begin(fastdiag.Rescale, "coefficient_loop", w.diagPar, fields)
	}
	for coefficient := 0; coefficient < ringQ.N(); coefficient++ {
		var timer fastdiag.Timer
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			timer = fastdiag.StartTimer()
		}
		var residues [MaxQPrefixWidth]uint64
		for row := 0; row < sourceRows; row++ {
			residues[row] = w.coeff.Coeffs[row][coefficient]
		}
		value := reconstructQPrefix(sourceRows, residues, w.q, w.inverse, w.modulus)
		magnitude, negative := centeredQPrefix(value, w.modulus[sourceRows-1], w.half[sourceRows-1])
		for step := 0; step < nbRescales; step++ {
			magnitude = roundedMagnitude192(magnitude, ringQ.SubRings[sourceLevel-step].Modulus)
			targetWidth, _ := QPrefixWidth(sourceLevel - step - 1)
			targetWidth = min(sourceRows, targetWidth)
			if cmp192(magnitude, w.half[targetWidth-1]) > 0 {
				if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
					reconstructNS += timer.ElapsedNS()
					loop.End(fastdiag.Fields{})
					fastdiag.Record(fastdiag.Rescale, "reconstruct_center_round_capacity", loop.Sequence(), fields, reconstructNS)
				}
				return &QPrefixCapacityError{Level: sourceLevel - step - 1, Component: component, Bound: uint192ToBigInt(magnitude), PrefixProduct: uint192ToBigInt(w.modulus[targetWidth-1])}
			}
		}
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			reconstructNS += timer.ElapsedNS()
			timer = fastdiag.StartTimer()
		}
		for row := 0; row < targetRows; row++ {
			staged.Coeffs[row][coefficient] = signedResidue192(magnitude, negative, w.q[row])
		}
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			materializeNS += timer.ElapsedNS()
		}
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		loop.End(fastdiag.Fields{})
		fastdiag.Record(fastdiag.Rescale, "reconstruct_center_round_capacity", loop.Sequence(), fields, reconstructNS)
		fastdiag.Record(fastdiag.Rescale, "residue_materialization", loop.Sequence(), fields, materializeNS)
	}
	return nil
}

func rescaleDiagFields(levelIn, levelOut, rowsIn, rowsOut int, component string, inPlace bool) fastdiag.Fields {
	return fastdiag.Fields{LevelIn: &levelIn, LevelOut: &levelOut, RowsIn: &rowsIn, RowsOut: &rowsOut, Component: component, InPlace: &inPlace}
}

func prefixToCoefficientRows(ringQ *ring.Ring, src ring.Poly, rows int, isNTT, isMontgomery bool, dst ring.Poly) error {
	if err := ValidatePrefixBacking(ringQ, rows, src, dst); err != nil {
		return err
	}
	for row := 0; row < rows; row++ {
		if isNTT {
			ringQ.SubRings[row].INTT(src.Coeffs[row], dst.Coeffs[row])
		} else {
			copy(dst.Coeffs[row], src.Coeffs[row])
		}
		if isMontgomery {
			ringQ.SubRings[row].IMForm(dst.Coeffs[row], dst.Coeffs[row])
		}
	}
	return nil
}

func reconstructQPrefix(rows int, residues [MaxQPrefixWidth]uint64, q, inverse [MaxQPrefixWidth]uint64, modulus [MaxQPrefixWidth]Uint192) Uint192 {
	switch rows {
	case 1:
		return Uint192{Lo: residues[0]}
	case 2:
		lo, hi := crtQ01(residues[0], residues[1], q[0], q[1], inverse[1])
		return Uint192{Lo: lo, Mid: hi}
	case 3:
		return crtQ012Prepared(residues[0], residues[1], residues[2], q[0], q[1], q[2], inverse[1], inverse[2], modulus[1])
	case 4:
		return crtQ0123Prepared(residues[0], residues[1], residues[2], residues[3], q[0], q[1], q[2], q[3], inverse[1], inverse[2], inverse[3], modulus[1], modulus[2])
	default:
		return Uint192{}
	}
}

func centeredQPrefix(value, modulus, half Uint192) (magnitude Uint192, negative bool) {
	if cmp192(value, half) > 0 {
		return sub192(modulus, value), true
	}
	return value, false
}

func uint192ToBigInt(value Uint192) *big.Int {
	out := new(big.Int).SetUint64(value.Hi)
	out.Lsh(out, 64).Add(out, new(big.Int).SetUint64(value.Mid))
	out.Lsh(out, 64).Add(out, new(big.Int).SetUint64(value.Lo))
	return out
}

func mul192By64(value Uint192, scalar uint64) (product Uint192, overflow bool) {
	hi0, lo0 := bits.Mul64(value.Lo, scalar)
	hi1, lo1 := bits.Mul64(value.Mid, scalar)
	hi2, lo2 := bits.Mul64(value.Hi, scalar)
	mid, carry := bits.Add64(hi0, lo1, 0)
	hi, carryOut := bits.Add64(hi1, lo2, carry)
	return Uint192{Lo: lo0, Mid: mid, Hi: hi}, hi2 != 0 || carryOut != 0
}

func add192(left, right Uint192) (sum Uint192, overflow bool) {
	lo, carry := bits.Add64(left.Lo, right.Lo, 0)
	mid, carry := bits.Add64(left.Mid, right.Mid, carry)
	hi, carry := bits.Add64(left.Hi, right.Hi, carry)
	return Uint192{Lo: lo, Mid: mid, Hi: hi}, carry != 0
}

func half192(value Uint192) Uint192 {
	return Uint192{Lo: (value.Lo >> 1) | (value.Mid << 63), Mid: (value.Mid >> 1) | (value.Hi << 63), Hi: value.Hi >> 1}
}

func sub192(left, right Uint192) Uint192 {
	lo, borrow := bits.Sub64(left.Lo, right.Lo, 0)
	mid, borrow := bits.Sub64(left.Mid, right.Mid, borrow)
	hi, _ := bits.Sub64(left.Hi, right.Hi, borrow)
	return Uint192{Lo: lo, Mid: mid, Hi: hi}
}

func cmp192(left, right Uint192) int {
	if left.Hi != right.Hi {
		if left.Hi < right.Hi {
			return -1
		}
		return 1
	}
	if left.Mid != right.Mid {
		if left.Mid < right.Mid {
			return -1
		}
		return 1
	}
	if left.Lo < right.Lo {
		return -1
	}
	if left.Lo > right.Lo {
		return 1
	}
	return 0
}

func mod192By64(value Uint192, modulus uint64) uint64 {
	remainder := value.Hi % modulus
	_, remainder = bits.Div64(remainder, value.Mid, modulus)
	_, remainder = bits.Div64(remainder, value.Lo, modulus)
	return remainder
}

func inverseMod(a, modulus uint64) (uint64, bool) {
	if modulus == 0 || a == 0 {
		return 0, false
	}
	oldR, r := int64(a), int64(modulus)
	oldT, t := int64(1), int64(0)
	for r != 0 {
		q := oldR / r
		oldR, r = r, oldR-q*r
		oldT, t = t, oldT-q*t
	}
	if oldR != 1 {
		return 0, false
	}
	if oldT < 0 {
		oldT += int64(modulus)
	}
	return uint64(oldT), true
}

func crtQ01(r0, r1, q0, q1, inverse uint64) (lo, hi uint64) {
	r0mod := r0 % q1
	var delta uint64
	if r1 >= r0mod {
		delta = r1 - r0mod
	} else {
		delta = q1 - (r0mod - r1)
	}
	prodHi, prodLo := bits.Mul64(delta, inverse)
	_, t := bits.Div64(prodHi, prodLo, q1)
	prodHi, prodLo = bits.Mul64(t, q0)
	lo, carry := bits.Add64(prodLo, r0, 0)
	hi, _ = bits.Add64(prodHi, 0, carry)
	return
}

func mod128By64(lo, hi, modulus uint64) uint64 {
	remainder := hi % modulus
	_, remainder = bits.Div64(remainder, lo, modulus)
	return remainder
}

func mul128By64(lo, hi, scalar uint64) Uint192 {
	hi0, lo0 := bits.Mul64(lo, scalar)
	hi1, lo1 := bits.Mul64(hi, scalar)
	mid, carry := bits.Add64(hi0, lo1, 0)
	hi2, _ := bits.Add64(hi1, 0, carry)
	return Uint192{Lo: lo0, Mid: mid, Hi: hi2}
}

func add128To192(lo, hi uint64, value Uint192) Uint192 {
	lo, carry := bits.Add64(value.Lo, lo, 0)
	mid, carry := bits.Add64(value.Mid, hi, carry)
	hi, _ = bits.Add64(value.Hi, 0, carry)
	return Uint192{Lo: lo, Mid: mid, Hi: hi}
}

func crtQ012Prepared(r0, r1, r2, q0, q1, q2, q0InvQ1, q01InvQ2 uint64, q01 Uint192) Uint192 {
	x01Lo, x01Hi := crtQ01(r0, r1, q0, q1, q0InvQ1)
	x01ModQ2 := mod128By64(x01Lo, x01Hi, q2)
	var delta uint64
	if r2 >= x01ModQ2 {
		delta = r2 - x01ModQ2
	} else {
		delta = q2 - (x01ModQ2 - r2)
	}
	prodHi, prodLo := bits.Mul64(delta, q01InvQ2)
	_, t := bits.Div64(prodHi, prodLo, q2)
	return add128To192(x01Lo, x01Hi, mul128By64(q01.Lo, q01.Mid, t))
}

func crtQ0123Prepared(r0, r1, r2, r3, q0, q1, q2, q3, q0InvQ1, q01InvQ2, q012InvQ3 uint64, q01, q012 Uint192) Uint192 {
	x012 := crtQ012Prepared(r0, r1, r2, q0, q1, q2, q0InvQ1, q01InvQ2, q01)
	x012ModQ3 := mod192By64(x012, q3)
	var delta uint64
	if r3 >= x012ModQ3 {
		delta = r3 - x012ModQ3
	} else {
		delta = q3 - (x012ModQ3 - r3)
	}
	prodHi, prodLo := bits.Mul64(delta, q012InvQ3)
	_, t := bits.Div64(prodHi, prodLo, q3)
	term, _ := mul192By64(q012, t)
	x, _ := add192(x012, term)
	return x
}

func roundedMagnitude192(value Uint192, divisor uint64) Uint192 {
	q2 := value.Hi / divisor
	remainder := value.Hi % divisor
	q1, remainder := bits.Div64(remainder, value.Mid, divisor)
	q0, remainder := bits.Div64(remainder, value.Lo, divisor)
	quotient := Uint192{Lo: q0, Mid: q1, Hi: q2}
	if remainder > divisor/2 {
		var carry uint64
		quotient.Lo, carry = bits.Add64(quotient.Lo, 1, 0)
		quotient.Mid, carry = bits.Add64(quotient.Mid, 0, carry)
		quotient.Hi, _ = bits.Add64(quotient.Hi, 0, carry)
	}
	return quotient
}

func signedResidue192(magnitude Uint192, negative bool, modulus uint64) uint64 {
	remainder := mod192By64(magnitude, modulus)
	if negative && remainder != 0 {
		return modulus - remainder
	}
	return remainder
}
