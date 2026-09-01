package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastN1ToN2 maps a ciphertext from a Standard ring of degree N1 to one of
// degree N2=2*N1. It is a pure q0/q1 ring-degree transformation: no
// evaluation key, key switch, CRT, or redistribution is involved.
//
// The input and output must be distinct ciphertexts with the same level and
// domain representation. q2 and above are deliberately left untouched in the
// output and are never read from the input.
func FastN1ToN2(ringN1, ringN2 *ring.Ring, ctIn, ctOut *rlwe.Ciphertext) error {
	return fastRingDegreeConversion(ringN1, ringN2, ctIn, ctOut, true)
}

// FastN2ToN1 maps a ciphertext from a Standard ring of degree N2=2*N1 to one
// of degree N1. In coefficient form it extracts every even coefficient. In
// NTT form it uses partial INTT, coefficient extraction, and partial NTT.
// q2 and above are deliberately left untouched in the output and are never
// read from the input.
func FastN2ToN1(ringN2, ringN1 *ring.Ring, ctIn, ctOut *rlwe.Ciphertext) error {
	return fastRingDegreeConversion(ringN2, ringN1, ctIn, ctOut, false)
}

func fastRingDegreeConversion(ringIn, ringOut *ring.Ring, ctIn, ctOut *rlwe.Ciphertext, expand bool) error {
	if ringIn == nil || ringOut == nil {
		return errors.New("input and output rings cannot be nil")
	}
	if ringIn.Type() != ring.Standard || ringOut.Type() != ring.Standard {
		return errors.New("Fast ring-degree conversion requires Standard rings")
	}
	if ringIn.Level() < 1 || ringOut.Level() < 1 {
		return errors.New("Fast ring-degree conversion requires q0 and q1")
	}
	if ringIn.Level() != ringOut.Level() {
		return errors.New("Fast ring-degree conversion requires equal ring levels")
	}
	if expand {
		if ringOut.N() != 2*ringIn.N() {
			return errors.New("FastN1ToN2 requires output degree twice the input degree")
		}
	} else if ringIn.N() != 2*ringOut.N() {
		return errors.New("FastN2ToN1 requires input degree twice the output degree")
	}
	for i := 0; i <= 1; i++ {
		if ringIn.SubRings[i].Modulus != ringOut.SubRings[i].Modulus {
			return errors.New("Fast ring-degree conversion requires matching q0/q1 moduli")
		}
	}
	if ctIn == nil || ctOut == nil {
		return errors.New("ctIn and ctOut cannot be nil")
	}
	if ctIn == ctOut {
		return errors.New("Fast ring-degree conversion requires distinct input and output ciphertexts")
	}
	if ctIn.MetaData == nil || ctOut.MetaData == nil {
		return errors.New("ctIn and ctOut metadata cannot be nil")
	}
	if ctIn.Degree() != 1 || ctOut.Degree() != 1 {
		return errors.New("Fast ring-degree conversion requires degree-one ciphertexts")
	}
	if ctIn.Level() != ctOut.Level() || ctIn.Level() != ringIn.Level() {
		return errors.New("Fast ring-degree conversion requires matching ciphertext and ring levels")
	}
	if ctIn.IsNTT != ctOut.IsNTT || ctIn.IsMontgomery != ctOut.IsMontgomery {
		return errors.New("Fast ring-degree conversion requires matching domains and Montgomery representations")
	}
	if ctIn.N() != ringIn.N() || ctOut.N() != ringOut.N() {
		return errors.New("ciphertext dimensions do not match ring-degree conversion")
	}
	for name, ct := range map[string]*rlwe.Ciphertext{"ctIn": ctIn, "ctOut": ctOut} {
		for d := 0; d <= 1; d++ {
			if len(ct.Value[d].Coeffs) < 2 || len(ct.Value[d].Coeffs[0]) != ct.N() || len(ct.Value[d].Coeffs[1]) != ct.N() {
				return fmt.Errorf("%s has invalid q0/q1 storage", name)
			}
		}
	}

	for d := 0; d <= 1; d++ {
		if ctIn.IsNTT {
			inCoeff := ring.NewPoly(ringIn.N(), 1)
			outCoeff := ring.NewPoly(ringOut.N(), 1)
			if err := FastPartialINTT(ringIn, ctIn.Value[d], inCoeff); err != nil {
				return fmt.Errorf("FastPartialINTT(%s): %w", componentName(d), err)
			}
			if err := mapRingDegreeCoefficient(inCoeff, outCoeff, expand); err != nil {
				return fmt.Errorf("map %s: %w", componentName(d), err)
			}
			if err := FastPartialNTT(ringOut, outCoeff, ctOut.Value[d]); err != nil {
				return fmt.Errorf("FastPartialNTT(%s): %w", componentName(d), err)
			}
		} else if err := mapRingDegreeCoefficient(ctIn.Value[d], ctOut.Value[d], expand); err != nil {
			return fmt.Errorf("map %s: %w", componentName(d), err)
		}
	}

	*ctOut.MetaData = *ctIn.MetaData
	ctOut.IsNTT = ctIn.IsNTT
	ctOut.IsMontgomery = ctIn.IsMontgomery
	return nil
}

func mapRingDegreeCoefficient(in, out ring.Poly, expand bool) error {
	if expand {
		for limb := 0; limb < 2; limb++ {
			for i, value := range in.Coeffs[limb] {
				out.Coeffs[limb][2*i] = value
				out.Coeffs[limb][2*i+1] = 0
			}
		}
	} else {
		for limb := 0; limb < 2; limb++ {
			for i := range out.Coeffs[limb] {
				out.Coeffs[limb][i] = in.Coeffs[limb][2*i]
			}
		}
	}
	return nil
}

func componentName(d int) string {
	if d == 0 {
		return "c0"
	}
	return "c1"
}
