package fast

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// Trace applies the ring-level Fast Trace using the legacy authoritative row
// count. It preserves ciphertext metadata and does not use evaluation keys or
// QP arithmetic. Boundaries that have explicitly materialized a wider
// Q-prefix should call TraceQPrefixRows instead.
func (eval *Evaluator) Trace(ctIn *rlwe.Ciphertext, logN int, opOut *rlwe.Ciphertext) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if ctIn == nil {
		return errors.New("ctIn and opOut cannot be nil")
	}
	rows := maintainedLimbCount(&eval.Parameters, ctIn.Level())
	return eval.traceQPrefixRows(ctIn, logN, rows, opOut, false)
}

// TraceQPrefixRows applies the ring-level Fast Trace to exactly rows explicit
// Q-prefix residues. It preserves the ciphertext scale and representation,
// supports in-place operation, and leaves all higher rows untouched.
func (eval *Evaluator) TraceQPrefixRows(ctIn *rlwe.Ciphertext, logN, rows int, opOut *rlwe.Ciphertext) error {
	return eval.traceQPrefixRows(ctIn, logN, rows, opOut, true)
}

func (eval *Evaluator) traceQPrefixRows(ctIn *rlwe.Ciphertext, logN, rows int, opOut *rlwe.Ciphertext, preserveHigherRows bool) error {
	if eval == nil {
		return errors.New("Fast evaluator cannot be nil")
	}
	if ctIn == nil || opOut == nil {
		return errors.New("ctIn and opOut cannot be nil")
	}
	if ctIn.MetaData == nil || opOut.MetaData == nil {
		return errors.New("ciphertext metadata cannot be nil")
	}
	if eval.Parameters.RingType() != ring.Standard {
		return fmt.Errorf("Fast Trace requires the Standard ring, got %s", eval.Parameters.RingType())
	}
	if logN < 0 || logN >= eval.Parameters.LogN() {
		return fmt.Errorf("Fast Trace logN must be in [0, %d), got %d", eval.Parameters.LogN(), logN)
	}
	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return errors.New("Fast Trace requires degree-one ciphertexts")
	}
	if ctIn.N() != eval.Parameters.N() || opOut.N() != eval.Parameters.N() {
		return errors.New("Fast Trace ciphertext dimensions do not match parameters")
	}
	if ctIn.Level() != opOut.Level() || ctIn.Level() < 1 {
		return errors.New("Fast Trace requires matching levels containing q0 and q1")
	}
	if ctIn.IsNTT != opOut.IsNTT || ctIn.IsMontgomery != opOut.IsMontgomery {
		return errors.New("Fast Trace requires matching NTT and Montgomery representations")
	}
	if !ctIn.IsNTT {
		return errors.New("Fast Trace requires NTT-domain ciphertexts")
	}
	for name, ct := range map[string]*rlwe.Ciphertext{"ctIn": ctIn, "opOut": opOut} {
		if len(ct.Value) != 2 {
			return fmt.Errorf("%s must have degree-one storage", name)
		}
	}
	if err := validatePrefixRows(eval.Parameters.RingQ(), ctIn.Level(), rows, ctIn.Value[0], ctIn.Value[1], opOut.Value[0], opOut.Value[1]); err != nil {
		return fmt.Errorf("Fast Trace Q-prefix storage: %w", err)
	}

	gap := 1 << (eval.Parameters.LogN() - logN - 1)
	if logN == 0 {
		gap <<= 1
	}
	var nInv *big.Int
	if gap > 1 {
		nInv = new(big.Int).SetUint64(uint64(gap))
		if nInv.ModInverse(nInv, eval.Parameters.RingQ().ModulusAtLevel[ctIn.Level()]) == nil {
			return fmt.Errorf("Fast Trace cannot invert gap %d modulo Q", gap)
		}
	}

	*opOut.MetaData = *ctIn.MetaData
	ringQ := eval.Parameters.RingQ().AtLevel(opOut.Level())
	if ctIn != opOut {
		for component := range ctIn.Value {
			copyPrefixRowsUnchecked(rows, ctIn.Value[component], opOut.Value[component])
		}
	}

	if gap <= 1 {
		return nil
	}
	var err error
	if preserveHigherRows {
		err = eval.MulIntegerQPrefixRows(opOut, nInv, rows, opOut)
	} else {
		err = eval.MulIntegerMaintained(opOut, nInv, opOut)
	}
	if err != nil {
		return fmt.Errorf("Fast Trace inverse normalization: %w", err)
	}

	apply := func(galEl uint64) error {
		for component := range opOut.Value {
			if err := eval.fastAutomorphismRows(ringQ, opOut.Value[component], eval.nttScratch[0], galEl, true, rows); err != nil {
				return err
			}
			for limb := 0; limb < rows; limb++ {
				ringQ.SubRings[limb].Add(opOut.Value[component].Coeffs[limb], eval.nttScratch[0].Coeffs[limb], opOut.Value[component].Coeffs[limb])
			}
		}
		return nil
	}

	for i := logN; i < eval.Parameters.LogN()-1; i++ {
		if err := apply(eval.Parameters.GaloisElement(1 << i)); err != nil {
			return fmt.Errorf("Fast Trace automorphism %d: %w", i, err)
		}
	}
	if logN == 0 {
		if err := apply(ringQ.NthRoot() - 1); err != nil {
			return fmt.Errorf("Fast Trace final order-two automorphism: %w", err)
		}
	}
	return nil
}
