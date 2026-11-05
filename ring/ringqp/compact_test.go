package ringqp

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func compactQPRing(t *testing.T) Ring {
	t.Helper()
	q, err := ring.NewRing(16, ring.Qi60[:6])
	require.NoError(t, err)
	p, err := ring.NewRing(16, ring.Pi60[:2])
	require.NoError(t, err)
	return Ring{q, p}
}

func compactQPInput(r Ring, qRows int, zeroP bool) Poly {
	p := r.NewPoly()
	for i, row := range p.Q.Coeffs {
		if i >= qRows {
			p.Q.Coeffs[i] = nil
			continue
		}
		for j := range row {
			row[j] = uint64(3 + i + j)
		}
	}
	if !zeroP {
		for i, row := range p.P.Coeffs {
			for j := range row {
				row[j] = uint64(9 + i + j)
			}
		}
	}
	return p
}

func TestCompactQPIndependentBasesAndModularOracle(t *testing.T) {
	r := compactQPRing(t)
	for _, zeroP := range []bool{false, true} {
		a, b := compactQPInput(r, 4, zeroP), compactQPInput(r, 6, zeroP)
		out := r.NewPoly()
		for repeat := 0; repeat < 3; repeat++ {
			r.MulCoeffsMontgomeryThenAdd(a, b, out)
			for side, params := range []*ring.Ring{r.RingQ, r.RingP} {
				x, y, z := a.Q, b.Q, out.Q
				if side == 1 {
					x, y, z = a.P, b.P, out.P
				}
				for i, row := range z.Coeffs {
					if side == 0 && i >= 4 {
						require.Empty(t, row)
						continue
					}
					require.Len(t, row, 16)
					q := new(big.Int).SetUint64(params.SubRings[i].Modulus)
					invR := new(big.Int).ModInverse(new(big.Int).Lsh(big.NewInt(1), 64), q)
					for j, v := range row {
						want := new(big.Int).Mul(new(big.Int).SetUint64(x.Coeffs[i][j]), new(big.Int).SetUint64(y.Coeffs[i][j]))
						want.Mul(want, invR)
						want.Mul(want, big.NewInt(int64(repeat+1)))
						want.Mod(want, q)
						require.Equal(t, want.Uint64(), v)
					}
				}
			}
		}
		// A full-source overwrite restores all Q rows, even on a compact receiver.
		r.Add(b, b, out)
		for side := 0; side < 2; side++ {
			x, z := b.Q, out.Q
			if side == 1 {
				x, z = b.P, out.P
			}
			for i, row := range z.Coeffs {
				require.Len(t, row, 16)
				for j, v := range row {
					require.Equal(t, 2*x.Coeffs[i][j], v)
				}
			}
		}
		index, err := ring.AutomorphismNTTIndex(16, r.RingQ.NthRoot(), 5)
		require.NoError(t, err)
		r.AutomorphismNTTWithIndex(a, index, out)
		for side := 0; side < 2; side++ {
			x, z := a.Q, out.Q
			if side == 1 {
				x, z = a.P, out.P
			}
			for i, row := range z.Coeffs {
				if side == 0 && i >= 4 {
					require.Empty(t, row)
					continue
				}
				for j, v := range row {
					require.Equal(t, x.Coeffs[i][index[j]], v)
				}
			}
		}
	}
}

func TestCompactQPPreflightBeforeEitherBaseWrites(t *testing.T) {
	r := compactQPRing(t)
	for _, operation := range []string{"add", "fma", "permutation"} {
		for _, kind := range []string{"P-hole", "P-short", "P-output-hole", "P-output-short", "P-alias", "cross-basis", "cross-headers", "index"} {
			t.Run(operation+"/"+kind, func(t *testing.T) {
				a, b, out := compactQPInput(r, 4, false), compactQPInput(r, 6, false), compactQPInput(r, 6, false)
				index, _ := ring.AutomorphismNTTIndex(16, r.RingQ.NthRoot(), 5)
				switch kind {
				case "P-hole":
					a.P.Coeffs[0] = nil
				case "P-short":
					a.P.Coeffs[1] = a.P.Coeffs[1][:15]
				case "P-output-hole":
					out.P.Coeffs[0] = nil
				case "P-output-short":
					out.P.Coeffs[1] = out.P.Coeffs[1][:15]
				case "P-alias":
					backing := make([]uint64, 17)
					a.P.Coeffs[0] = backing[:16]
					out.P.Coeffs[0] = backing[1:]
				case "cross-basis":
					out.Q.Coeffs[0] = a.P.Coeffs[0]
				case "cross-headers":
					out.P.Coeffs = out.Q.Coeffs[:2]
				case "index":
					if operation != "permutation" {
						t.Skip("only a permutation has an index")
					}
					index = out.P.Coeffs[0]
					for i := range index {
						index[i] = uint64(i)
					}
				}
				// Raw snapshot also covers deliberately malformed output shapes.
				before := &Poly{Q: ring.Poly{Coeffs: out.Q.Coeffs.CopyNew()}, P: ring.Poly{Coeffs: out.P.Coeffs.CopyNew()}}
				require.Panics(t, func() {
					switch operation {
					case "add":
						r.Add(a, b, out)
					case "fma":
						r.MulCoeffsMontgomeryThenAdd(a, b, out)
					case "permutation":
						r.AutomorphismNTTWithIndex(a, index, out)
					}
				})
				require.True(t, out.Equal(before))
			})
		}
	}
}

func TestCompactQPLegalInPlaceAndRepeatedZero(t *testing.T) {
	r := compactQPRing(t)
	a, b := compactQPInput(r, 4, false), compactQPInput(r, 6, false)
	before := a.CopyNew()
	r.Add(a, b, a)
	for i, row := range a.Q.Coeffs[:4] {
		for j, v := range row {
			require.Equal(t, before.Q.Coeffs[i][j]+b.Q.Coeffs[i][j], v)
		}
	}
	for i := 4; i < 6; i++ {
		require.Empty(t, a.Q.Coeffs[i])
	}
	for i, row := range a.P.Coeffs {
		for j, v := range row {
			require.Equal(t, before.P.Coeffs[i][j]+b.P.Coeffs[i][j], v)
		}
	}
	for repeat := 0; repeat < 3; repeat++ {
		a.Q.Zero()
		a.P.Zero()
		r.MulCoeffsMontgomeryThenAdd(b, b, a)
		for i := 4; i < 6; i++ {
			require.Empty(t, a.Q.Coeffs[i])
		}
	}
}
