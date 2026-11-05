package ring

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
)

func TestZeroPModDownUsesExistingConstantsAndModularOracle(t *testing.T) {
	q := compactTestRing(t)
	p, err := NewRing(q.N(), Pi60[:3])
	require.NoError(t, err)
	be := NewBasisExtender(q, p)
	for _, rows := range []int{1, 4, 6} {
		for _, levelP := range []int{0, 2} {
			in := compactFixture(16, 5, rows)
			for _, mont := range []bool{false, true} {
				source := *in.CopyNew()
				if mont {
					q.MForm(source, source)
				}
				before := source.CopyNew()
				for _, alias := range []bool{false, true} {
					out := q.NewPoly()
					if alias {
						out = *source.CopyNew()
					}
					input := source
					if alias {
						input = out
					}
					require.NoError(t, be.ModDownQPtoQNTTZeroP(5, levelP, rows, input, p.AtLevel(levelP).NewPoly(), out))
					for i, row := range out.Coeffs {
						if i >= rows {
							require.Empty(t, row)
							continue
						}
						modulus := new(big.Int).SetUint64(q.SubRings[i].Modulus)
						inverse := new(big.Int).ModInverse(p.ModulusAtLevel[levelP], modulus)
						for j, value := range row {
							want := new(big.Int).Mul(new(big.Int).SetUint64(before.Coeffs[i][j]), inverse)
							want.Mod(want, modulus)
							require.Equal(t, want.Uint64(), value)
						}
					}
				}
				require.True(t, source.Equal(before))
			}
		}
	}
}

func TestZeroPModDownRejectsBeforeWrites(t *testing.T) {
	q := compactTestRing(t)
	p, err := NewRing(16, Pi60[:2])
	require.NoError(t, err)
	be := NewBasisExtender(q, p)
	for _, kind := range []string{"nonzero-P", "missing-P", "missing-Q", "zero-width", "short-output", "partial-alias", "cross-P-alias"} {
		t.Run(kind, func(t *testing.T) {
			in, ps, out := compactFixture(16, 5, 4), p.NewPoly(), q.NewPoly()
			switch kind {
			case "nonzero-P":
				ps.Coeffs[0][0] = 1
			case "missing-P":
				ps.Coeffs[1] = nil
			case "missing-Q":
				in.Coeffs[3] = nil
			case "zero-width":
				in = compactFixture(16, 5, 0)
			case "short-output":
				out.Coeffs[4] = out.Coeffs[4][:15]
			case "partial-alias":
				x := make([]uint64, 17)
				in.Coeffs[0] = x[:16]
				out.Coeffs[0] = x[1:]
			case "cross-P-alias":
				out.Coeffs[0] = ps.Coeffs[0]
			}
			before := out.Coeffs.CopyNew()
			require.Error(t, be.ModDownQPtoQNTTZeroP(5, 1, 4, in, ps, out))
			require.Equal(t, before, out.Coeffs)
		})
	}
}
