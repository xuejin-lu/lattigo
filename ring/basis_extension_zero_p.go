package ring

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring/internal/compact"
)

// ModDownQPtoQNTTZeroP is the exact zero-P specialization of ModDownQPtoQNTT.
// It checks P itself, uses the existing Montgomery P-inverse constants, and
// preserves the logical Q level while retaining only the requested Q rows.
// It does not reconstruct an integer or assert any CKKS capacity condition.
func (be *BasisExtender) ModDownQPtoQNTTZeroP(levelQ, levelP, rows int, p1Q, p1P, p2Q Poly) error {
	if be == nil || be.ringQ == nil || be.ringP == nil ||
		levelQ < 0 || levelQ > be.ringQ.MaxLevel() || levelP < 0 || levelP > be.ringP.MaxLevel() ||
		rows < 1 || rows > levelQ+1 {
		return fmt.Errorf("zero-P ModDown: invalid rings, levels or Q row count")
	}
	n := be.ringQ.N()
	if be.ringP.N() != n {
		return fmt.Errorf("zero-P ModDown: Q/P dimensions differ")
	}
	width, err := compact.ValidateWidth(p1Q.Coeffs, levelQ, n)
	if err != nil {
		return err
	}
	if width < rows {
		return fmt.Errorf("zero-P ModDown: missing required Q rows")
	}
	if _, err = compact.ValidateWidth(p2Q.Coeffs, levelQ, n); err != nil {
		return err
	}
	if rows < levelQ+1 && p2Q.Level() != levelQ {
		return fmt.Errorf("zero-P ModDown: compact output level must match")
	}
	if width, err = compact.ValidateWidth(p1P.Coeffs, levelP, n); err != nil {
		return err
	}
	if width != levelP+1 {
		return fmt.Errorf("zero-P ModDown: missing required P rows")
	}
	for _, row := range p1P.Coeffs[:levelP+1] {
		for _, value := range row {
			if value != 0 {
				return fmt.Errorf("zero-P ModDown: P input is nonzero")
			}
		}
	}
	if err = compact.ValidateAliases(p2Q.Coeffs, levelQ, true, false, p1Q.Coeffs); err != nil {
		return err
	}
	for _, dst := range p2Q.Coeffs {
		for _, src := range p1P.Coeffs {
			if compact.Overlap(dst, src) {
				return fmt.Errorf("zero-P ModDown: output aliases P input")
			}
		}
	}
	// No writes or pool use precede complete validation. Exact in-place Q is legal.
	compact.Prepare(p2Q.Coeffs, levelQ, n, rows)
	for row, subring := range be.ringQ.SubRings[:rows] {
		subring.MulScalarMontgomery(p1Q.Coeffs[row], be.modDownConstantsPtoQ[levelP][row], p2Q.Coeffs[row])
	}
	return nil
}
