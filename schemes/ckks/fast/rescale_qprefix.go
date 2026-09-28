package fast

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/internal/fastdiag"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// rescaleNQPrefix implements Standard's sequential Rescale rounding over the
// authoritative Q-prefix. It preflights every exact result before mutating the
// destination so capacity failures are transactional, including in-place use.
func (eval *Evaluator) rescaleNQPrefix(op0 *rlwe.Ciphertext, nbRescales, sourceRows int, opOut *rlwe.Ciphertext) error {
	var wholeSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		wholeSpan = fastdiag.Begin(fastdiag.Rescale, "rescale", 0, fastdiag.Input(op0, sourceRows).Merge(fastdiag.InPlace(op0 == opOut)).Merge(fastdiag.RepetitionCount(nbRescales)))
	}
	if eval == nil || op0 == nil || opOut == nil {
		return errors.New("Fast Rescale evaluator and ciphertexts cannot be nil")
	}
	if op0.MetaData == nil || opOut.MetaData == nil {
		return errors.New("Fast Rescale metadata cannot be nil")
	}
	if nbRescales <= 0 || op0.Level() < nbRescales {
		return errors.New("invalid number of Fast rescale levels")
	}
	if op0.N() != eval.Parameters.N() || opOut.N() != eval.Parameters.N() {
		return errors.New("ciphertext dimensions do not match Fast evaluator parameters")
	}
	if err := validateFastRescaleDomain(op0); err != nil {
		return err
	}
	ringQ := eval.Parameters.RingQ()
	if ringQ == nil || op0.Level() > ringQ.Level() {
		return fmt.Errorf("Fast Rescale input Level %d exceeds configured ring Level", op0.Level())
	}
	maxSourceRows, err := QPrefixWidth(op0.Level())
	if err != nil {
		return err
	}
	if sourceRows < 1 || sourceRows > maxSourceRows {
		return fmt.Errorf("explicit Fast Rescale source width %d is outside [1,%d] at Level %d", sourceRows, maxSourceRows, op0.Level())
	}
	targetLevel := op0.Level() - nbRescales
	maxTargetRows, err := QPrefixWidth(targetLevel)
	if err != nil {
		return err
	}
	targetRows := min(sourceRows, maxTargetRows)

	scratch := &eval.rescaleScratch
	if scratch.coeff.N() != ringQ.N() || len(scratch.coeff.Coeffs) < sourceRows {
		eval.rescaleScratch = newFastRescaleScratch(ringQ)
		scratch = &eval.rescaleScratch
	}
	if err := scratch.validateWidth(sourceRows); err != nil {
		return err
	}
	if sourceRows > len(scratch.q) || sourceRows > len(scratch.coeff.Coeffs) || targetRows > len(scratch.result.Coeffs) {
		return fmt.Errorf("Fast Rescale scratch does not cover source/target widths %d/%d", sourceRows, targetRows)
	}
	if len(op0.Value) == 0 {
		return errors.New("Fast Rescale ciphertext has no components")
	}
	for component := range op0.Value {
		if err := validatePrefixRows(ringQ, op0.Level(), sourceRows, op0.Value[component]); err != nil {
			return fmt.Errorf("Fast Rescale component %d: %w", component, err)
		}
	}
	for step := 0; step < nbRescales; step++ {
		logicalLevel := op0.Level() - step
		if logicalLevel >= len(ringQ.SubRings) || ringQ.SubRings[logicalLevel] == nil || ringQ.SubRings[logicalLevel].Modulus < 2 {
			return fmt.Errorf("Fast Rescale logical divisor q%d is unavailable", logicalLevel)
		}
	}

	// Pass one computes only into evaluator-owned scratch and rejects any
	// intermediate prefix overflow before opOut or an in-place op0 is touched.
	var preflightSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		preflightSpan = fastdiag.Begin(fastdiag.Rescale, "preflight", wholeSpan.Sequence(), fastdiag.Input(op0, sourceRows).Merge(fastdiag.InPlace(op0 == opOut)).Merge(fastdiag.RepetitionCount(len(op0.Value))))
		scratch.diagParent = preflightSpan.Sequence()
	}
	for component := range op0.Value {
		if err := eval.rescaleQPrefixComponent(op0.Value[component], nil, op0.Level(), nbRescales, sourceRows, targetRows, component, op0.IsMontgomery, op0 == opOut, scratch); err != nil {
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
		scratch.diagParent = materializationSpan.Sequence()
	}

	if opOut != op0 {
		Resize(opOut, op0.Degree(), targetLevel, eval.Parameters.N())
	}
	for component := range op0.Value {
		if err := eval.rescaleQPrefixComponent(op0.Value[component], &opOut.Value[component], op0.Level(), nbRescales, sourceRows, targetRows, component, op0.IsMontgomery, op0 == opOut, scratch); err != nil {
			return err
		}
	}
	*opOut.MetaData = *op0.MetaData
	opOut.Scale = op0.Scale
	for step := 0; step < nbRescales; step++ {
		opOut.Scale = opOut.Scale.Div(rlwe.NewScale(ringQ.SubRings[op0.Level()-step].Modulus))
	}
	if opOut == op0 {
		Resize(opOut, op0.Degree(), targetLevel, eval.Parameters.N())
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		rows, _ := QPrefixWidth(opOut.Level())
		materializationSpan.End(fastdiag.Output(opOut, rows))
		wholeSpan.End(fastdiag.Output(opOut, rows))
	}
	return nil
}

func (eval *Evaluator) rescaleQPrefixComponent(src ring.Poly, dst *ring.Poly, sourceLevel, nbRescales, sourceRows, targetRows, component int, montgomery, inPlace bool, scratch *fastRescaleScratch) error {
	ringQ := eval.Parameters.RingQ()
	var fields fastdiag.Fields
	var phaseSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		componentName := fmt.Sprintf("c%d", component)
		fields = fastDiagRescaleFields(sourceLevel, sourceLevel-nbRescales, sourceRows, targetRows, componentName, inPlace)
		phaseSpan = fastdiag.Begin(fastdiag.Rescale, "prefix_to_coefficient", scratch.diagParent, fields)
	}
	if err := prefixToCoefficientRows(ringQ, src, sourceRows, true, montgomery, scratch.coeff); err != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			phaseSpan.End(fastdiag.Fields{})
		}
		return err
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		phaseSpan.End(fastdiag.Fields{})
	}

	var reconstructNS, residueMaterializationNS int64
	var coefficientLoopSpan fastdiag.Span
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		coefficientLoopSpan = fastdiag.Begin(fastdiag.Rescale, "coefficient_loop", scratch.diagParent, fields)
	}
	for coefficient := 0; coefficient < ringQ.N(); coefficient++ {
		var segment fastdiag.Timer
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			segment = fastdiag.StartTimer()
		}
		var residues [MaxQPrefixWidth]uint64
		for row := 0; row < sourceRows; row++ {
			residues[row] = scratch.coeff.Coeffs[row][coefficient]
		}
		value := reconstructQPrefix(sourceRows, residues, scratch)
		magnitude, negative := centeredQPrefix(value, scratch.modulus[sourceRows-1], scratch.half[sourceRows-1])
		for step := 0; step < nbRescales; step++ {
			magnitude = roundedMagnitude192(magnitude, ringQ.SubRings[sourceLevel-step].Modulus)
			targetWidth, _ := QPrefixWidth(sourceLevel - step - 1)
			targetWidth = min(sourceRows, targetWidth)
			if cmp192(magnitude, scratch.half[targetWidth-1]) > 0 {
				if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
					reconstructNS += segment.ElapsedNS()
					coefficientLoopSpan.End(fastdiag.Fields{})
					fastdiag.Record(fastdiag.Rescale, "reconstruct_center_round_capacity", coefficientLoopSpan.Sequence(), fields, reconstructNS)
				}
				return &QPrefixCapacityError{
					Level:         sourceLevel - step - 1,
					Component:     component,
					Bound:         uint192ToBigInt(magnitude),
					PrefixProduct: uint192ToBigInt(scratch.modulus[targetWidth-1]),
				}
			}
		}
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			reconstructNS += segment.ElapsedNS()
		}
		if dst != nil {
			if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
				segment = fastdiag.StartTimer()
			}
			for row := 0; row < targetRows; row++ {
				scratch.result.Coeffs[row][coefficient] = signedResidue192(magnitude, negative, scratch.q[row])
			}
			if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
				residueMaterializationNS += segment.ElapsedNS()
			}
		}
	}
	if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
		coefficientLoopSpan.End(fastdiag.Fields{})
		fastdiag.Record(fastdiag.Rescale, "reconstruct_center_round_capacity", coefficientLoopSpan.Sequence(), fields, reconstructNS)
		if dst != nil {
			fastdiag.Record(fastdiag.Rescale, "residue_materialization", coefficientLoopSpan.Sequence(), fields, residueMaterializationNS)
		}
	}

	if dst != nil {
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			phaseSpan = fastdiag.Begin(fastdiag.Rescale, "ntt_montgomery_restore", scratch.diagParent, fields)
		}
		for row := 0; row < targetRows; row++ {
			ringQ.SubRings[row].NTT(scratch.result.Coeffs[row], dst.Coeffs[row])
			if montgomery {
				ringQ.SubRings[row].MForm(dst.Coeffs[row], dst.Coeffs[row])
			}
		}
		if fastdiag.Enabled && fastdiag.Selected(fastdiag.Rescale) {
			phaseSpan.End(fastdiag.Fields{})
		}
	}
	return nil
}

func fastDiagRescaleFields(levelIn, levelOut, rowsIn, rowsOut int, component string, inPlace bool) fastdiag.Fields {
	return fastdiag.Fields{
		LevelIn: &levelIn, LevelOut: &levelOut, RowsIn: &rowsIn, RowsOut: &rowsOut,
		Component: component, InPlace: &inPlace,
	}
}

func prefixToCoefficientRows(ringQ *ring.Ring, src ring.Poly, rows int, isNTT, isMontgomery bool, dst ring.Poly) error {
	if err := validatePrefixBacking(ringQ, rows, src, dst); err != nil {
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

func reconstructQPrefix(rows int, residues [MaxQPrefixWidth]uint64, scratch *fastRescaleScratch) uint192 {
	switch rows {
	case 1:
		return uint192{lo: residues[0]}
	case 2:
		lo, hi := crtQ01(residues[0], residues[1], scratch.q[0], scratch.q[1], scratch.inverse[1])
		return uint192{lo: lo, mid: hi}
	case 3:
		return crtQ012(residues[0], residues[1], residues[2], scratch.q[0], scratch.q[1], scratch.q[2], scratch.inverse[1], scratch.inverse[2])
	case 4:
		return crtQ0123(residues[0], residues[1], residues[2], residues[3], scratch.q[0], scratch.q[1], scratch.q[2], scratch.q[3], scratch.inverse[1], scratch.inverse[2], scratch.inverse[3])
	default:
		return uint192{}
	}
}

func centeredQPrefix(value, modulus, half uint192) (magnitude uint192, negative bool) {
	if cmp192(value, half) > 0 {
		return sub192(modulus, value), true
	}
	return value, false
}

func uint192ToBigInt(value uint192) *big.Int {
	out := new(big.Int).SetUint64(value.hi)
	out.Lsh(out, 64)
	out.Add(out, new(big.Int).SetUint64(value.mid))
	out.Lsh(out, 64)
	out.Add(out, new(big.Int).SetUint64(value.lo))
	return out
}
