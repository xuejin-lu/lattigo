package ring

import "github.com/tuneinsight/lattigo/v6/ring/internal/compact"

// prepareCompactOperation preserves availability without selecting any CKKS
// width. All validation precedes row-header changes and existing SubRing calls.
// Shallow copies sharing row headers are in-place aliases. Unsupported partial
// overlap is rejected. Rows above the Ring's requested level remain untouched.
func (r Ring) prepareCompactOperation(out Poly, accumulate, exactAlias bool, inputs ...Poly) int {
	coeffs := make([][][]uint64, len(inputs))
	for i := range inputs {
		coeffs[i] = inputs[i].Coeffs
	}
	rows := compact.Plan(r.N(), r.level, out.Coeffs, accumulate, exactAlias, coeffs...)
	compact.Prepare(out.Coeffs, r.level, r.N(), rows)
	return rows
}

func (r Ring) prepareCompactTransform(in, out Poly) int {
	rows := compact.PlanTransform(r.N(), r.level, out.Coeffs, in.Coeffs)
	compact.Prepare(out.Coeffs, r.level, r.N(), rows)
	return rows
}
