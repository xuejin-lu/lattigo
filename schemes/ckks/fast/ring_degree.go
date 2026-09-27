package fast

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// FastN1ToN2 maps a ciphertext from a Standard ring of degree N1 to one of
// degree N2=2*N1 using the legacy maintained-row authority. It is a pure
// ring-degree transformation: no evaluation key, key switch, CRT, or
// redistribution is involved.
func FastN1ToN2(ringN1, ringN2 *ring.Ring, ctIn, ctOut *rlwe.Ciphertext) error {
	return FastN1ToN2QPrefixRows(ringN1, ringN2, ctIn, ctOut, maintainedLimbCountForRing(ringN1))
}

// FastN2ToN1 maps a ciphertext from a Standard ring of degree N2=2*N1 to one
// of degree N1 using the legacy maintained-row authority. In coefficient form
// it extracts every even coefficient. In NTT form it uses partial INTT,
// coefficient extraction, and partial NTT.
func FastN2ToN1(ringN2, ringN1 *ring.Ring, ctIn, ctOut *rlwe.Ciphertext) error {
	return FastN2ToN1QPrefixRows(ringN2, ringN1, ctIn, ctOut, maintainedLimbCountForRing(ringN2))
}

// FastN1ToN2QPrefixRows maps exactly rows authoritative Q-prefix rows from N1
// to N2. Rows outside the requested prefix are never read or written.
func FastN1ToN2QPrefixRows(ringN1, ringN2 *ring.Ring, ctIn, ctOut *rlwe.Ciphertext, rows int) error {
	return fastRingDegreeConversionRows(ringN1, ringN2, ctIn, ctOut, true, rows)
}

// FastN2ToN1QPrefixRows maps exactly rows authoritative Q-prefix rows from N2
// to N1. Rows outside the requested prefix are never read or written.
func FastN2ToN1QPrefixRows(ringN2, ringN1 *ring.Ring, ctIn, ctOut *rlwe.Ciphertext, rows int) error {
	return fastRingDegreeConversionRows(ringN2, ringN1, ctIn, ctOut, false, rows)
}

func fastRingDegreeConversionRows(ringIn, ringOut *ring.Ring, ctIn, ctOut *rlwe.Ciphertext, expand bool, rows int) error {
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
	maxRows := min(MaxQPrefixWidth, ringIn.Level()+1)
	if rows < 1 || rows > maxRows {
		return fmt.Errorf("Fast ring-degree conversion rows %d must be in [1,%d]", rows, maxRows)
	}
	if expand {
		if ringOut.N() != 2*ringIn.N() {
			return errors.New("FastN1ToN2 requires output degree twice the input degree")
		}
	} else if ringIn.N() != 2*ringOut.N() {
		return errors.New("FastN2ToN1 requires input degree twice the output degree")
	}
	for i := 0; i < rows; i++ {
		if ringIn.SubRings[i].Modulus != ringOut.SubRings[i].Modulus {
			return fmt.Errorf("Fast ring-degree conversion requires matching q%d moduli", i)
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
	width, err := QPrefixWidth(ctIn.Level())
	if err != nil || rows > width {
		return fmt.Errorf("Fast ring-degree conversion rows %d exceed the Q-prefix width at level %d", rows, ctIn.Level())
	}
	if ctIn.IsNTT != ctOut.IsNTT || ctIn.IsMontgomery != ctOut.IsMontgomery {
		return errors.New("Fast ring-degree conversion requires matching domains and Montgomery representations")
	}
	if ctIn.N() != ringIn.N() || ctOut.N() != ringOut.N() {
		return errors.New("ciphertext dimensions do not match ring-degree conversion")
	}
	for name, ct := range map[string]*rlwe.Ciphertext{"ctIn": ctIn, "ctOut": ctOut} {
		for d := 0; d <= 1; d++ {
			if len(ct.Value[d].Coeffs) < rows {
				return fmt.Errorf("%s has invalid maintained storage", name)
			}
			for limb := 0; limb < rows; limb++ {
				if len(ct.Value[d].Coeffs[limb]) != ct.N() {
					return fmt.Errorf("%s has invalid maintained storage", name)
				}
			}
		}
	}

	for d := 0; d <= 1; d++ {
		if ctIn.IsNTT {
			inCoeff := ring.NewPoly(ringIn.N(), rows-1)
			outCoeff := ring.NewPoly(ringOut.N(), rows-1)
			if err := FastPartialINTTRows(ringIn, ctIn.Value[d], inCoeff, rows); err != nil {
				return fmt.Errorf("FastPartialINTT(%s): %w", componentName(d), err)
			}
			if err := mapRingDegreeCoefficientRows(inCoeff, outCoeff, expand, rows); err != nil {
				return fmt.Errorf("map %s: %w", componentName(d), err)
			}
			if err := FastPartialNTTRows(ringOut, outCoeff, ctOut.Value[d], rows); err != nil {
				return fmt.Errorf("FastPartialNTT(%s): %w", componentName(d), err)
			}
		} else if err := mapRingDegreeCoefficientRows(ctIn.Value[d], ctOut.Value[d], expand, rows); err != nil {
			return fmt.Errorf("map %s: %w", componentName(d), err)
		}
	}

	*ctOut.MetaData = *ctIn.MetaData
	ctOut.IsNTT = ctIn.IsNTT
	ctOut.IsMontgomery = ctIn.IsMontgomery
	return nil
}

func mapRingDegreeCoefficientRows(in, out ring.Poly, expand bool, rows int) error {
	if rows < 1 || len(in.Coeffs) < rows || len(out.Coeffs) < rows {
		return errors.New("ring-degree coefficient mapping has insufficient rows")
	}
	if expand {
		for limb := 0; limb < rows; limb++ {
			if len(out.Coeffs[limb]) != 2*len(in.Coeffs[limb]) {
				return fmt.Errorf("ring-degree expansion q%d has mismatched coefficient lengths", limb)
			}
			for i, value := range in.Coeffs[limb] {
				out.Coeffs[limb][2*i] = value
				out.Coeffs[limb][2*i+1] = 0
			}
		}
	} else {
		for limb := 0; limb < rows; limb++ {
			if len(in.Coeffs[limb]) != 2*len(out.Coeffs[limb]) {
				return fmt.Errorf("ring-degree contraction q%d has mismatched coefficient lengths", limb)
			}
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
