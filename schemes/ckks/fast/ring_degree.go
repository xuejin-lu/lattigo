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
	if ringIn.Level() < 0 || ringOut.Level() < 0 {
		return errors.New("Fast ring-degree conversion requires q0")
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
	maintained := ringIn.Level() + 1
	if maintained > 2 {
		maintained = 2
	}
	for i := 0; i < maintained; i++ {
		if ringIn.SubRings[i].Modulus != ringOut.SubRings[i].Modulus {
			return errors.New("Fast ring-degree conversion requires matching maintained moduli")
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
			if len(ct.Value[d].Coeffs) < maintained || len(ct.Value[d].Coeffs[0]) != ct.N() {
				return fmt.Errorf("%s has invalid maintained storage", name)
			}
			for limb := 1; limb < maintained; limb++ {
				if len(ct.Value[d].Coeffs[limb]) != ct.N() {
					return fmt.Errorf("%s has invalid maintained storage", name)
				}
			}
		}
	}

	for d := 0; d <= 1; d++ {
		if ctIn.IsNTT {
			inCoeff := ring.NewPoly(ringIn.N(), maintained)
			outCoeff := ring.NewPoly(ringOut.N(), maintained)
			if maintained == 1 {
				ringIn.SubRings[0].INTT(ctIn.Value[d].Coeffs[0], inCoeff.Coeffs[0])
			} else if err := FastPartialINTT(ringIn, ctIn.Value[d], inCoeff); err != nil {
				return fmt.Errorf("FastPartialINTT(%s): %w", componentName(d), err)
			}
			if err := mapRingDegreeCoefficient(inCoeff, outCoeff, expand); err != nil {
				return fmt.Errorf("map %s: %w", componentName(d), err)
			}
			if maintained == 1 {
				ringOut.SubRings[0].NTT(outCoeff.Coeffs[0], ctOut.Value[d].Coeffs[0])
			} else if err := FastPartialNTT(ringOut, outCoeff, ctOut.Value[d]); err != nil {
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
		maintained := len(in.Coeffs)
		if maintained > 2 {
			maintained = 2
		}
		for limb := 0; limb < maintained; limb++ {
			for i, value := range in.Coeffs[limb] {
				out.Coeffs[limb][2*i] = value
				out.Coeffs[limb][2*i+1] = 0
			}
		}
	} else {
		maintained := len(out.Coeffs)
		if maintained > 2 {
			maintained = 2
		}
		for limb := 0; limb < maintained; limb++ {
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
