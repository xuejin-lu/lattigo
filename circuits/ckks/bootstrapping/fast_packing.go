package bootstrapping

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// ensureFastPackingTables lazily creates the Standard-ring monomial tables
// used by the packing boundary. The tables contain only q0/q1 rows; the
// exponent schedule is the same as the Standard evaluator's tables.
func (eval *FastEvaluator) ensureFastPackingTables() error {
	if eval == nil {
		return errors.New("Fast Bootstrap evaluator cannot be nil")
	}
	if eval.fastPackingInitialized {
		return eval.fastPackingErr
	}
	eval.fastPackingInitialized = true

	paramsN1 := eval.Parameters.ResidualParameters
	paramsN2 := eval.Parameters.BootstrappingParameters
	if paramsN2.RingQ() == nil {
		eval.fastPackingErr = errors.New("Fast packing requires BootstrappingParameters")
		return eval.fastPackingErr
	}
	if paramsN1.RingQ() == nil && paramsN1.N() != paramsN2.N() {
		eval.fastPackingErr = errors.New("Fast packing requires ResidualParameters for a ring-degree boundary")
		return eval.fastPackingErr
	}
	if paramsN2.RingType() != ring.Standard || (paramsN1.RingQ() != nil && paramsN1.RingType() != ring.Standard) {
		eval.fastPackingErr = errors.New("Fast packing requires Standard rings")
		return eval.fastPackingErr
	}
	if paramsN2.N() < paramsN1.N() {
		eval.fastPackingErr = fmt.Errorf("Fast packing does not support N2 < N1: N1=%d N2=%d", paramsN1.N(), paramsN2.N())
		return eval.fastPackingErr
	}
	if paramsN2.N() > 2*paramsN1.N() {
		eval.fastPackingErr = fmt.Errorf("Fast packing supports at most N2=2*N1: N1=%d N2=%d", paramsN1.N(), paramsN2.N())
		return eval.fastPackingErr
	}
	if paramsN2.RingQ().MaxLevel() < 0 {
		eval.fastPackingErr = errors.New("Fast packing requires q0 in BootstrappingParameters")
		return eval.fastPackingErr
	}
	if paramsN2.N() != paramsN1.N() && paramsN1.RingQ().MaxLevel() < 0 {
		eval.fastPackingErr = errors.New("Fast packing requires q0 in ResidualParameters")
		return eval.fastPackingErr
	}
	if paramsN1.RingQ() != nil {
		maintained := utils.Min(paramsN1.RingQ().MaxLevel(), paramsN2.RingQ().MaxLevel()) + 1
		maintained = utils.Min(maintained, 2)
		for limb := 0; limb < maintained; limb++ {
			if paramsN1.RingQ().SubRings[limb].Modulus != paramsN2.RingQ().SubRings[limb].Modulus {
				eval.fastPackingErr = errors.New("Fast packing requires matching maintained moduli across N1/N2")
				return eval.fastPackingErr
			}
		}
	}

	paramsN2Q01 := paramsN2.RingQ().AtLevel(utils.Min(1, paramsN2.RingQ().MaxLevel()))
	eval.xPow2N2 = rlwe.GenXPow2NTT(paramsN2Q01, paramsN2.LogN(), false)
	eval.xPow2InvN2 = rlwe.GenXPow2NTT(paramsN2Q01, paramsN2.LogN(), true)

	if paramsN1.N() != paramsN2.N() {
		paramsN1Q01 := paramsN1.RingQ().AtLevel(utils.Min(1, paramsN1.RingQ().MaxLevel()))
		// Standard uses paramsN2.LogN() for the forward N1 table and
		// paramsN1.LogN() for its inverse table.
		eval.xPow2N1 = rlwe.GenXPow2NTT(paramsN1Q01, paramsN2.LogN(), false)
		eval.xPow2InvN1 = rlwe.GenXPow2NTT(paramsN1Q01, paramsN1.LogN(), true)
	}
	return nil
}

func validateFastPackingCiphertext(ct *rlwe.Ciphertext, params ckks.Parameters) error {
	if ct == nil || ct.MetaData == nil {
		return errors.New("Fast packing ciphertext and metadata cannot be nil")
	}
	if params.RingQ() == nil || params.RingType() != ring.Standard {
		return errors.New("Fast packing requires a Standard parameter ring")
	}
	if ct.N() != params.N() {
		return fmt.Errorf("Fast packing ciphertext degree %d does not match parameters degree %d", ct.N(), params.N())
	}
	if ct.Degree() != 1 {
		return fmt.Errorf("Fast packing requires degree-one ciphertexts, got degree %d", ct.Degree())
	}
	if ct.Level() < 0 || ct.Level() > params.RingQ().MaxLevel() {
		return fmt.Errorf("Fast packing requires ciphertext level in [0, %d], got %d", params.RingQ().MaxLevel(), ct.Level())
	}
	if !ct.IsNTT {
		return errors.New("Fast packing requires NTT-domain ciphertexts")
	}
	maintained := maintainedLimbs(ct.Level())
	for d := 0; d <= 1; d++ {
		if len(ct.Value) <= d || len(ct.Value[d].Coeffs) < maintained ||
			len(ct.Value[d].Coeffs[0]) != params.N() {
			return errors.New("Fast packing requires maintained ciphertext storage")
		}
		for limb := 1; limb < maintained; limb++ {
			if len(ct.Value[d].Coeffs[limb]) != params.N() {
				return errors.New("Fast packing requires maintained ciphertext storage")
			}
		}
	}
	return nil
}

func validateFastPackingSlice(cts []rlwe.Ciphertext, params ckks.Parameters) error {
	if len(cts) == 0 {
		return errors.New("Fast packing requires a non-empty ciphertext slice")
	}
	if err := validateFastPackingCiphertext(&cts[0], params); err != nil {
		return err
	}
	level := cts[0].Level()
	logSlots := cts[0].LogSlots()
	for i := 1; i < len(cts); i++ {
		if err := validateFastPackingCiphertext(&cts[i], params); err != nil {
			return fmt.Errorf("ciphertext %d: %w", i, err)
		}
		if cts[i].Level() != level {
			return fmt.Errorf("ciphertext %d level %d does not match level %d", i, cts[i].Level(), level)
		}
		if cts[i].IsMontgomery != cts[0].IsMontgomery {
			return fmt.Errorf("ciphertext %d Montgomery representation does not match ciphertext 0", i)
		}
		if cts[i].IsBatched != cts[0].IsBatched || cts[i].IsBitReversed != cts[0].IsBitReversed {
			return fmt.Errorf("ciphertext %d plaintext representation does not match ciphertext 0", i)
		}
		if cts[i].LogSlots() != logSlots {
			return fmt.Errorf("ciphertext %d LogSlots %d does not match LogSlots %d", i, cts[i].LogSlots(), logSlots)
		}
	}
	return nil
}

func copyFastPackingCiphertext(params ckks.Parameters, src *rlwe.Ciphertext) rlwe.Ciphertext {
	dst := fastckks.NewCiphertext(params, 1, src.Level())
	*dst.MetaData = *src.MetaData
	maintained := maintainedLimbs(src.Level())
	for d := 0; d <= 1; d++ {
		for limb := 0; limb < maintained; limb++ {
			copy(dst.Value[d].Coeffs[limb], src.Value[d].Coeffs[limb])
		}
	}
	return *dst
}

func (eval *FastEvaluator) fastPack(cts []rlwe.Ciphertext, ctxt packingContext, xPow2 []ring.Poly) ([]rlwe.Ciphertext, error) {
	if ctxt.Params == nil {
		return nil, errors.New("Fast packing context parameters cannot be nil")
	}
	if ctxt.LogSlots < 0 || ctxt.LogSlots > ctxt.LogMaxDimensions.Cols || ctxt.NbPackedCTs < 1 {
		return nil, errors.New("invalid Fast packing context")
	}
	if err := validateFastPackingSlice(cts, *ctxt.Params); err != nil {
		return nil, err
	}
	if len(xPow2) < ctxt.LogMaxDimensions.Cols-ctxt.LogSlots {
		return nil, errors.New("Fast packing monomial table is too short")
	}

	packed := make([]rlwe.Ciphertext, len(cts))
	for i := range cts {
		packed[i] = copyFastPackingCiphertext(*ctxt.Params, &cts[i])
	}
	logPackCTs := ctxt.LogMaxDimensions.Cols - ctxt.LogSlots
	logGap := ctxt.Params.LogMaxSlots() - ctxt.LogSlots - 1
	maintained := maintainedLimbs(packed[0].Level())
	for i := 0; i < logPackCTs; i++ {
		for j := 0; j < len(packed)>>1; j++ {
			even := &packed[j*2]
			odd := &packed[j*2+1]
			if even.Level() != odd.Level() {
				return nil, errors.New("Fast packing requires equal levels within each pair")
			}
			monomialIndex := logGap - i
			if monomialIndex < 0 || monomialIndex >= len(xPow2) {
				return nil, fmt.Errorf("Fast packing monomial index %d is out of range", monomialIndex)
			}
			monomial := xPow2[monomialIndex]
			ringQ := ctxt.Params.RingQ()
			for d := 0; d <= 1; d++ {
				for limb := 0; limb < maintained; limb++ {
					ringQ.SubRings[limb].MulCoeffsMontgomeryThenAdd(
						odd.Value[d].Coeffs[limb], monomial.Coeffs[limb], even.Value[d].Coeffs[limb])
				}
			}
			packed[j] = *even
		}
		if len(packed)&1 == 1 {
			packed[len(packed)>>1] = packed[len(packed)-1]
		}
		if len(packed)&1 == 0 {
			packed = packed[:len(packed)>>1]
		} else {
			packed = packed[:len(packed)>>1+1]
		}
	}
	for i := range packed {
		packed[i].LogDimensions = ctxt.LogMaxDimensions
	}
	return packed, nil
}

func (eval *FastEvaluator) fastUnpack(ct *rlwe.Ciphertext, ctxt packingContext, xPow2Inv []ring.Poly) ([]rlwe.Ciphertext, error) {
	if ctxt.Params == nil {
		return nil, errors.New("Fast unpacking context parameters cannot be nil")
	}
	if ctxt.LogSlots < 0 || ctxt.LogSlots > ctxt.LogMaxDimensions.Cols || ctxt.NbPackedCTs < 1 {
		return nil, errors.New("invalid Fast unpacking context")
	}
	if err := validateFastPackingCiphertext(ct, *ctxt.Params); err != nil {
		return nil, err
	}
	logPackCTs := ctxt.LogMaxDimensions.Cols - ctxt.LogSlots
	if len(xPow2Inv) < logPackCTs {
		return nil, errors.New("Fast unpacking monomial table is too short")
	}
	n := utils.Min(ctxt.NbPackedCTs, 1<<logPackCTs)
	if n < 1 {
		return nil, errors.New("Fast unpacking context has no ciphertexts to unpack")
	}
	cts := make([]rlwe.Ciphertext, n)
	for i := range cts {
		cts[i] = copyFastPackingCiphertext(*ctxt.Params, ct)
	}

	logGap := ctxt.Params.LogMaxSlots() - ctxt.LogSlots - 1
	maintained := maintainedLimbs(ct.Level())
	for i := 0; i < utils.Min(bitsLen64(uint64(n-1)), logPackCTs); i++ {
		step := 1 << (i + 1)
		monomialIndex := logGap - i
		if monomialIndex < 0 || monomialIndex >= len(xPow2Inv) {
			return nil, fmt.Errorf("Fast unpacking monomial index %d is out of range", monomialIndex)
		}
		monomial := xPow2Inv[monomialIndex]
		for j := 0; j < n; j += step {
			for k := step >> 1; k < step && j+k < n; k++ {
				for d := 0; d <= 1; d++ {
					for limb := 0; limb < maintained; limb++ {
						ctPoly := cts[j+k].Value[d].Coeffs[limb]
						ctxt.Params.RingQ().SubRings[limb].MulCoeffsMontgomery(ctPoly, monomial.Coeffs[limb], ctPoly)
					}
				}
			}
		}
	}
	return cts, nil
}

func maintainedLimbs(level int) int {
	return utils.Min(level+1, 2)
}

func bitsLen64(value uint64) int {
	var length int
	for value != 0 {
		length++
		value >>= 1
	}
	return length
}

// PackAndSwitchN1ToN2 packs Fast ciphertexts in N1, maps them to N2 when
// necessary, and packs them again in N2. It does not apply evaluation keys.
func (eval *FastEvaluator) PackAndSwitchN1ToN2(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, *packingContext, *packingContext, error) {
	if err := eval.ensureFastPackingTables(); err != nil {
		return nil, nil, nil, err
	}
	paramsN1 := eval.Parameters.ResidualParameters
	paramsN2 := eval.Parameters.BootstrappingParameters
	if paramsN1.RingQ() == nil {
		paramsN1 = paramsN2
	}
	if err := validateFastPackingSlice(cts, paramsN1); err != nil {
		return nil, nil, nil, err
	}

	var packN1 *packingContext
	if paramsN1.N() != paramsN2.N() {
		packN1 = &packingContext{Params: &eval.Parameters.ResidualParameters, LogMaxDimensions: paramsN1.LogMaxDimensions(), LogSlots: cts[0].LogSlots(), NbPackedCTs: len(cts)}
		if eval.Parameters.LogMaxSlots() < paramsN1.LogMaxSlots() {
			packN1.LogMaxDimensions = eval.Parameters.LogMaxDimensions()
		}
		var err error
		if cts, err = eval.fastPack(cts, *packN1, eval.xPow2N1); err != nil {
			return nil, nil, nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN1: %w", err)
		}
		for i := range cts {
			out := fastckks.NewCiphertext(paramsN2, 1, cts[i].Level())
			out.IsNTT = cts[i].IsNTT
			out.IsMontgomery = cts[i].IsMontgomery
			if err := fastckks.FastN1ToN2(paramsN1.RingQ().AtLevel(cts[i].Level()), paramsN2.RingQ().AtLevel(cts[i].Level()), &cts[i], out); err != nil {
				return nil, nil, nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: FastN1ToN2: %w", err)
			}
			cts[i] = *out
		}
	}

	packN2 := &packingContext{Params: &eval.Parameters.BootstrappingParameters, LogMaxDimensions: eval.Parameters.LogMaxDimensions(), LogSlots: cts[0].LogSlots(), NbPackedCTs: len(cts)}
	var err error
	if cts, err = eval.fastPack(cts, *packN2, eval.xPow2N2); err != nil {
		return nil, nil, nil, fmt.Errorf("cannot PackAndSwitchN1ToN2: PackN2: %w", err)
	}
	return cts, packN1, packN2, nil
}

// UnpackAndSwitchN2ToN1 reverses PackAndSwitchN1ToN2 using only maintained
// q0/q1 residues and the pure Fast ring-degree mapping.
func (eval *FastEvaluator) UnpackAndSwitchN2ToN1(cts []rlwe.Ciphertext, ctxtN1, ctxtN2 *packingContext) ([]rlwe.Ciphertext, error) {
	if err := eval.ensureFastPackingTables(); err != nil {
		return nil, err
	}
	if len(cts) == 0 {
		return nil, errors.New("Fast unpacking requires a non-empty ciphertext slice")
	}
	if ctxtN2 == nil {
		return nil, errors.New("Fast unpacking requires an N2 packing context")
	}
	var ctsOut []rlwe.Ciphertext
	logSlots := ctxtN2.LogSlots
	for i := range cts {
		unpacked, err := eval.fastUnpack(&cts[i], *ctxtN2, eval.xPow2InvN2)
		if err != nil {
			return nil, fmt.Errorf("cannot UnpackAndSwitchN2ToN1: UnpackN2: %w", err)
		}
		ctsOut = append(ctsOut, unpacked...)
		ctxtN2.NbPackedCTs -= len(unpacked)
	}

	if ctxtN1 != nil {
		if ctxtN1.Params == nil || ctxtN1.NbPackedCTs < 1 {
			return nil, errors.New("invalid N1 packing context")
		}
		logSlots = ctxtN1.LogSlots
		ctsN1 := make([]rlwe.Ciphertext, 0, len(ctsOut))
		paramsN1 := eval.Parameters.ResidualParameters
		paramsN2 := eval.Parameters.BootstrappingParameters
		for i := range ctsOut {
			out := fastckks.NewCiphertext(paramsN1, 1, ctsOut[i].Level())
			out.IsNTT = ctsOut[i].IsNTT
			out.IsMontgomery = ctsOut[i].IsMontgomery
			if err := fastckks.FastN2ToN1(paramsN2.RingQ().AtLevel(ctsOut[i].Level()), paramsN1.RingQ().AtLevel(ctsOut[i].Level()), &ctsOut[i], out); err != nil {
				return nil, fmt.Errorf("cannot UnpackAndSwitchN2ToN1: FastN2ToN1: %w", err)
			}
			unpacked, err := eval.fastUnpack(out, *ctxtN1, eval.xPow2InvN1)
			if err != nil {
				return nil, fmt.Errorf("cannot UnpackAndSwitchN2ToN1: UnpackN1: %w", err)
			}
			ctsN1 = append(ctsN1, unpacked...)
			ctxtN1.NbPackedCTs -= len(unpacked)
		}
		ctsOut = ctsN1
	}
	for i := range ctsOut {
		ctsOut[i].LogDimensions.Cols = logSlots
	}
	return ctsOut, nil
}
