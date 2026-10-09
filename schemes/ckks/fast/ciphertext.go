package fast

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

// NewCiphertext allocates a Fast-local ciphertext shape. The logical level is
// retained in the coefficient-row headers, while only the Q-prefix rows own
// N-sized backing arrays. Rows outside the prefix are nil and must not be
// passed to Standard or full-RNS APIs.
func NewCiphertext(params rlwe.ParameterProvider, degree, level int) *rlwe.Ciphertext {
	return fastcore.NewCompactCiphertext(params, degree, level)
}

// Resize changes Fast ciphertext logical level/degree without materializing
// dormant rows. It is the compact counterpart of rlwe.Element.Resize.
func Resize(ct *rlwe.Ciphertext, degree, level, n int) {
	fastcore.ResizeCompactCiphertext(ct, degree, level, n)
}
