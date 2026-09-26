package fast

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

// NewCiphertext allocates a Fast-local ciphertext shape. The logical level is
// retained in the coefficient-row headers, while only the Q-prefix rows own
// N-sized backing arrays. Rows outside the prefix are nil and must not be
// passed to Standard or full-RNS APIs.
func NewCiphertext(params rlwe.ParameterProvider, degree, level int) *rlwe.Ciphertext {
	n := params.GetRLWEParameters().N()
	prefixWidth := qPrefixWidthOrPanic(level)
	polys := make([]ring.Poly, degree+1)
	for i := range polys {
		polys[i] = compactPoly(n, level, prefixWidth)
	}
	ct, err := rlwe.NewCiphertextAtLevelFromPoly(level, polys)
	if err != nil {
		panic(err)
	}
	ct.IsNTT = params.GetRLWEParameters().NTTFlag()
	return ct
}

// Resize changes Fast ciphertext logical level/degree without materializing
// dormant rows. It is the compact counterpart of rlwe.Element.Resize.
func Resize(ct *rlwe.Ciphertext, degree, level, n int) {
	if ct == nil || level < 0 || degree < 0 {
		panic("invalid Fast ciphertext resize")
	}

	prefixWidth := qPrefixWidthOrPanic(level)
	for i := range ct.Value {
		coeffs := ct.Value[i].Coeffs
		if level < len(coeffs)-1 {
			ct.Value[i].Coeffs = coeffs[:level+1]
		} else if level > len(coeffs)-1 {
			grown := make([][]uint64, level+1)
			copy(grown, coeffs)
			ct.Value[i].Coeffs = grown
		}
		ensureQPrefixRows(ct.Value[i], level, n, prefixWidth)
	}

	if degree < len(ct.Value)-1 {
		ct.Value = ct.Value[:degree+1]
	} else {
		for len(ct.Value) < degree+1 {
			ct.Value = append(ct.Value, compactPoly(n, level, prefixWidth))
		}
	}
}

func compactPoly(n, level, prefixWidth int) ring.Poly {
	coeffs := make([][]uint64, level+1)
	for limb := 0; limb < len(coeffs) && limb < prefixWidth; limb++ {
		coeffs[limb] = make([]uint64, n)
	}
	return ring.Poly{Coeffs: coeffs}
}

func ensureQPrefixRows(poly ring.Poly, level, n, prefixWidth int) {
	for limb := range poly.Coeffs {
		if limb < prefixWidth {
			if len(poly.Coeffs[limb]) != n {
				poly.Coeffs[limb] = make([]uint64, n)
			}
		} else {
			poly.Coeffs[limb] = nil
		}
	}
}
