package ckks

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

// fastCKKSZeroSecretCiphertextLevel performs the structural checks needed
// before RotateNew allocates an output from the input's logical level.
func fastCKKSZeroSecretCiphertextLevel(ct *rlwe.Ciphertext) (int, error) {
	if ct == nil {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: ciphertext cannot be nil")
	}
	if len(ct.Value) != 2 {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: expected degree-one ciphertext, got %d polynomial components", len(ct.Value))
	}
	if len(ct.Value[0].Coeffs) == 0 || len(ct.Value[1].Coeffs) == 0 {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: ciphertext has no active Q rows")
	}
	level := len(ct.Value[0].Coeffs) - 1
	if len(ct.Value[1].Coeffs) != level+1 {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: c1 logical row count does not match c0")
	}
	return level, nil
}

func inspectFastCKKSZeroSecretRotateCiphertext(name string, ct *rlwe.Ciphertext, n, maxLevel int) (int, error) {
	level, err := fastCKKSZeroSecretCiphertextLevel(ct)
	if err != nil {
		return 0, fmt.Errorf("%s: %w", name, err)
	}
	if ct.MetaData == nil {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: %s metadata cannot be nil", name)
	}
	if level > maxLevel {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: %s level %d exceeds configured maximum %d", name, level, maxLevel)
	}
	if ct.IsMontgomery && !ct.IsNTT {
		return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: %s Montgomery representation requires NTT domain", name)
	}
	for component := range ct.Value {
		poly := ct.Value[component]
		if len(poly.Coeffs) != level+1 {
			return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: %s component %d logical rows do not match level %d", name, component, level)
		}
		for row := 0; row <= level; row++ {
			if len(poly.Coeffs[row]) != n {
				return 0, fmt.Errorf("cannot Fast CKKS zero-secret Rotate: %s component %d q%d row is not fully materialized at N=%d", name, component, row, n)
			}
		}
	}
	return level, nil
}

// rotateFastCKKSZeroSecret applies the shared ring-level Fast permutation to
// the complete active-Q ciphertext. It deliberately bypasses Galois-key lookup
// only for the Fast zero-secret capability and rejects unsupported inputs
// before the shared core can write any output row.
func (eval Evaluator) rotateFastCKKSZeroSecret(ctIn *rlwe.Ciphertext, galEl uint64, ctOut *rlwe.Ciphertext) error {
	if ctIn == nil || ctOut == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret Rotate: input and output ciphertexts cannot be nil")
	}
	params := eval.GetParameters()
	if params.RingType() != ring.Standard {
		return fmt.Errorf("cannot Fast CKKS zero-secret Rotate: requires the Standard ring, got %s", params.RingType())
	}
	if eval.automorphismCore == nil {
		return fmt.Errorf("cannot Fast CKKS zero-secret Rotate: shared automorphism core is unavailable")
	}
	level, err := inspectFastCKKSZeroSecretRotateCiphertext("input", ctIn, params.N(), params.MaxLevel())
	if err != nil {
		return err
	}
	outLevel, err := inspectFastCKKSZeroSecretRotateCiphertext("output", ctOut, params.N(), params.MaxLevel())
	if err != nil {
		return err
	}
	if level != outLevel {
		return fmt.Errorf("cannot Fast CKKS zero-secret Rotate: input level %d and output level %d must match", level, outLevel)
	}
	if ctIn.IsNTT != ctOut.IsNTT || ctIn.IsMontgomery != ctOut.IsMontgomery {
		return fmt.Errorf("cannot Fast CKKS zero-secret Rotate: input and output domains/representations must match")
	}
	for row := 0; row <= level; row++ {
		modulus := params.RingQ().SubRings[row].Modulus
		for coefficient, value := range ctIn.Value[0].Coeffs[row] {
			if value >= modulus {
				return fmt.Errorf("cannot Fast CKKS zero-secret Rotate: input c0 q%d coefficient %d is not a canonical residue", row, coefficient)
			}
		}
		for coefficient, value := range ctIn.Value[1].Coeffs[row] {
			if value != 0 {
				return fmt.Errorf("cannot Fast CKKS zero-secret Rotate: input c1 q%d coefficient %d is nonzero", row, coefficient)
			}
		}
	}

	rows := level + 1
	ringQ := params.RingQ().AtLevel(level)
	if err := eval.automorphismCore.ApplyRows(ringQ, galEl, ctIn.IsNTT, rows,
		fastcore.PolynomialPair{Input: ctIn.Value[0], Output: ctOut.Value[0]},
		fastcore.PolynomialPair{Input: ctIn.Value[1], Output: ctOut.Value[1]},
	); err != nil {
		return fmt.Errorf("cannot apply shared Fast automorphism core: %w", err)
	}
	if !ctIn.IsNTT {
		// The reused coefficient permutation represents -0 as q, matching the
		// ring primitive. Canonicalize that equivalent residue in the public
		// zero-secret contract so c1 remains exactly zero and rows stay in [0,q).
		for row := 0; row <= level; row++ {
			modulus := ringQ.SubRings[row].Modulus
			for component := range ctOut.Value {
				for coefficient, value := range ctOut.Value[component].Coeffs[row] {
					if value == modulus {
						ctOut.Value[component].Coeffs[row][coefficient] = 0
					}
				}
			}
		}
	}
	*ctOut.MetaData = *ctIn.MetaData
	return nil
}
