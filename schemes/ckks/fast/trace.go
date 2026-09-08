package fast

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// Trace applies the ring-level Fast Trace to a degree-one ciphertext. It
// preserves the ciphertext scale and representation, reads and writes only
// q0/q1, and does not use evaluation keys or QP arithmetic.
func (eval *Evaluator) Trace(ctIn *rlwe.Ciphertext, logN int, opOut *rlwe.Ciphertext) error {
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
		for component := range ct.Value {
			if len(ct.Value[component].Coeffs) < 2 ||
				len(ct.Value[component].Coeffs[0]) != eval.Parameters.N() ||
				len(ct.Value[component].Coeffs[1]) != eval.Parameters.N() {
				return fmt.Errorf("%s has invalid q0/q1 storage", name)
			}
		}
	}

	*opOut.MetaData = *ctIn.MetaData
	if ctIn != opOut {
		for component := range ctIn.Value {
			copy(opOut.Value[component].Coeffs[0], ctIn.Value[component].Coeffs[0])
			copy(opOut.Value[component].Coeffs[1], ctIn.Value[component].Coeffs[1])
		}
	}

	gap := 1 << (eval.Parameters.LogN() - logN - 1)
	if logN == 0 {
		gap <<= 1
	}
	if gap <= 1 {
		return nil
	}

	ringQ := eval.Parameters.RingQ().AtLevel(opOut.Level())
	nInv := new(big.Int).SetUint64(uint64(gap))
	if nInv.ModInverse(nInv, ringQ.ModulusAtLevel[opOut.Level()]) == nil {
		return fmt.Errorf("Fast Trace cannot invert gap %d modulo Q", gap)
	}
	if err := eval.MulIntegerMaintained(opOut, nInv, opOut); err != nil {
		return fmt.Errorf("Fast Trace inverse normalization: %w", err)
	}

	apply := func(galEl uint64) error {
		for component := range opOut.Value {
			if err := eval.fastAutomorphism(ringQ, opOut.Value[component], eval.nttScratch[0], galEl, true); err != nil {
				return err
			}
			for limb := 0; limb < 2; limb++ {
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
