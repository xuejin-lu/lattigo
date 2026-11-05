package ringqp

import "github.com/tuneinsight/lattigo/v6/ring/internal/compact"

// preflightCompact checks both bases before either underlying Ring can write.
// Availability is independent in Q/P; this does not assert zero P or a lift.
func (r Ring) preflightCompact(out Poly, accumulate, exactAlias bool, inputs ...Poly) {
	if r.RingQ != nil && r.RingP != nil && r.RingQ.N() != r.RingP.N() {
		panic("ringqp compact: Q/P ring dimensions differ")
	}
	var qInputs, pInputs [][][]uint64
	for _, in := range inputs {
		qInputs = append(qInputs, in.Q.Coeffs)
		pInputs = append(pInputs, in.P.Coeffs)
	}
	if r.RingQ != nil {
		compact.Plan(r.RingQ.N(), r.RingQ.Level(), out.Q.Coeffs, accumulate, exactAlias, qInputs...)
	}
	if r.RingP != nil {
		compact.Plan(r.RingP.N(), r.RingP.Level(), out.P.Coeffs, accumulate, exactAlias, pInputs...)
	}
	// Also prevent Q commits from changing an input or output of the P phase.
	compact.Disjoint(out.Q.Coeffs, out.P.Coeffs)
	for _, in := range inputs {
		compact.Disjoint(out.Q.Coeffs, in.P.Coeffs)
		compact.Disjoint(out.P.Coeffs, in.Q.Coeffs)
	}
}
