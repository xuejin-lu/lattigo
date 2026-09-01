package fast

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastMulQ01 multiplies two coefficient-domain polynomials using only q0 and
// q1 for the NTT-domain computation. The result is materialized in all Q
// limbs in coefficient domain. The bound is an explicit precondition: every
// coefficient of the exact negacyclic product must be in [-bound, bound].
//
// The input dormant limbs are never read or modified. The output dormant
// limbs are populated only by reducing the reconstructed product coefficients
// modulo the corresponding qi.
func FastMulQ01(ringQ *ring.Ring, p1, p2, p3 ring.Poly, bound *big.Int) error {
	if err := validateFastMul(ringQ, p1, p2, p3, bound); err != nil {
		return err
	}

	level := ringQ.Level()
	cCoeff, err := fastMulQ01Core(ringQ, p1, p2)
	if err != nil {
		return err
	}

	backend := NewRNSBackend(ringQ, level)
	coeffs, err := backend.ReconstructQ0Q1(&cCoeff)
	if err != nil {
		return fmt.Errorf("ReconstructQ0Q1(product): %w", err)
	}
	for i, coeff := range coeffs {
		if err := backend.CheckCoefficientBounds(coeff, bound); err != nil {
			return fmt.Errorf("product coefficient %d exceeds bound: %w", i, err)
		}
	}

	if err := backend.Redistribute(coeffs, &p3); err != nil {
		return fmt.Errorf("Redistribute(product): %w", err)
	}
	return nil
}

// FastMulQ01Authoritative multiplies two coefficient-domain polynomials using
// only q0 and q1. It is the production Fast contract: q0 and q1 are written,
// while q2 and above are deliberately neither read nor written. No CRT,
// bound check or redistribution is performed.
func FastMulQ01Authoritative(ringQ *ring.Ring, p1, p2, p3 ring.Poly) error {
	if err := validateFastMulAuthoritative(ringQ, p1, p2, p3); err != nil {
		return err
	}
	product, err := fastMulQ01Core(ringQ.AtLevel(ringQ.Level()), p1, p2)
	if err != nil {
		return err
	}
	copy(p3.Coeffs[0], product.Coeffs[0])
	copy(p3.Coeffs[1], product.Coeffs[1])
	return nil
}

func fastMulQ01Core(ringQ *ring.Ring, p1, p2 ring.Poly) (ring.Poly, error) {
	aNTT := ring.NewPoly(ringQ.N(), 1)
	bNTT := ring.NewPoly(ringQ.N(), 1)
	cNTT := ring.NewPoly(ringQ.N(), 1)
	cCoeff := ring.NewPoly(ringQ.N(), 1)

	if err := FastPartialNTT(ringQ, p1, aNTT); err != nil {
		return ring.Poly{}, fmt.Errorf("FastPartialNTT(p1): %w", err)
	}
	if err := FastPartialNTT(ringQ, p2, bNTT); err != nil {
		return ring.Poly{}, fmt.Errorf("FastPartialNTT(p2): %w", err)
	}

	// MulCoeffsMontgomery expects one operand in Montgomery form. This is the
	// same representation transition used by the Standard CKKS multiplication.
	ringQ.SubRings[0].MForm(aNTT.Coeffs[0], aNTT.Coeffs[0])
	ringQ.SubRings[1].MForm(aNTT.Coeffs[1], aNTT.Coeffs[1])
	ringQ.SubRings[0].MulCoeffsMontgomery(aNTT.Coeffs[0], bNTT.Coeffs[0], cNTT.Coeffs[0])
	ringQ.SubRings[1].MulCoeffsMontgomery(aNTT.Coeffs[1], bNTT.Coeffs[1], cNTT.Coeffs[1])

	if err := FastPartialINTT(ringQ, cNTT, cCoeff); err != nil {
		return ring.Poly{}, fmt.Errorf("FastPartialINTT(product): %w", err)
	}
	return cCoeff, nil
}

func validateFastMulAuthoritative(ringQ *ring.Ring, p1, p2, p3 ring.Poly) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Level() < 1 {
		return errors.New("FastMulQ01Authoritative requires q0 and q1")
	}
	for name, p := range map[string]ring.Poly{"p1": p1, "p2": p2} {
		if p.N() != ringQ.N() || p.Level() < ringQ.Level() {
			return fmt.Errorf("%s dimensions or level are insufficient", name)
		}
		if len(p.Coeffs) < 2 || len(p.Coeffs[0]) != ringQ.N() || len(p.Coeffs[1]) != ringQ.N() {
			return fmt.Errorf("%s q0/q1 coefficient slices are invalid", name)
		}
	}
	if p3.N() != ringQ.N() || p3.Level() < 1 {
		return errors.New("p3 dimensions or level are insufficient")
	}
	if len(p3.Coeffs) < 2 || len(p3.Coeffs[0]) != ringQ.N() || len(p3.Coeffs[1]) != ringQ.N() {
		return errors.New("p3 q0/q1 coefficient slices are invalid")
	}
	return nil
}

func validateFastMul(ringQ *ring.Ring, p1, p2, p3 ring.Poly, bound *big.Int) error {
	if ringQ == nil {
		return errors.New("ringQ cannot be nil")
	}
	if ringQ.Level() < 1 {
		return errors.New("FastMulQ01 requires q0 and q1")
	}
	if bound == nil {
		return errors.New("coefficient bound cannot be nil")
	}
	if bound.Sign() < 0 {
		return errors.New("coefficient bound cannot be negative")
	}
	if p1.Level() != ringQ.Level() || p2.Level() != ringQ.Level() || p3.Level() != ringQ.Level() {
		return fmt.Errorf("input and output levels must equal ringQ.Level()=%d", ringQ.Level())
	}
	if p1.N() != ringQ.N() || p2.N() != ringQ.N() || p3.N() != ringQ.N() {
		return fmt.Errorf("input and output dimensions must match ringQ.N()=%d", ringQ.N())
	}
	for name, p := range map[string]ring.Poly{"p1": p1, "p2": p2, "p3": p3} {
		for i := 0; i <= ringQ.Level(); i++ {
			if len(p.Coeffs[i]) != ringQ.N() {
				return fmt.Errorf("%s coefficient slice %d must have length %d", name, i, ringQ.N())
			}
		}
	}
	return nil
}
