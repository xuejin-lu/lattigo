package fast

import (
	"fmt"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

type prefixTestCase struct {
	level int
	rows  int
}

func prefixTestCases() []prefixTestCase {
	return []prefixTestCase{{0, 1}, {1, 2}, {2, 3}, {3, 1}, {3, 2}, {3, 3}, {3, 4}}
}

func fillPrefixPoly(poly ring.Poly, ringQ *ring.Ring, salt uint64) {
	for row := range poly.Coeffs {
		q := ringQ.SubRings[row].Modulus
		for i := range poly.Coeffs[row] {
			poly.Coeffs[row][i] = (salt + uint64(97*row+13*i+1)) % q
		}
	}
}

func fillPrefixCiphertext(ct *rlwe.Ciphertext, ringQ *ring.Ring, salt uint64) {
	for component := range ct.Value {
		fillPrefixPoly(ct.Value[component], ringQ, salt+uint64(component*1009))
	}
}

func fillPrefixSentinel(poly ring.Poly, ringQ *ring.Ring, salt uint64) {
	for row := range poly.Coeffs {
		q := ringQ.SubRings[row].Modulus
		for i := range poly.Coeffs[row] {
			poly.Coeffs[row][i] = (q - 1 - (salt+uint64(17*row+i))%q) % q
		}
	}
}

func clonePrefixPoly(src ring.Poly) ring.Poly {
	dst := ring.NewPoly(src.N(), src.Level())
	for row := range src.Coeffs {
		copy(dst.Coeffs[row], src.Coeffs[row])
	}
	return dst
}

func requirePrefixRowsEqual(t *testing.T, want, got ring.Poly, rows int) {
	t.Helper()
	for row := 0; row < rows; row++ {
		require.Equal(t, want.Coeffs[row], got.Coeffs[row], "q%d", row)
	}
}

func prefixMulMod(x, y, modulus uint64) uint64 {
	product := new(big.Int).Mul(new(big.Int).SetUint64(x), new(big.Int).SetUint64(y))
	return product.Mod(product, new(big.Int).SetUint64(modulus)).Uint64()
}

func prefixAddMod(x, y, modulus uint64) uint64 {
	sum := new(big.Int).Add(new(big.Int).SetUint64(x), new(big.Int).SetUint64(y))
	return sum.Mod(sum, new(big.Int).SetUint64(modulus)).Uint64()
}

func TestPrefixCompleteAddSubCopyNegateAndIntegerMultiply(t *testing.T) {
	params := testFastCKKSParameters(t)
	ringQ := params.RingQ()
	eval := NewEvaluator(params)

	for _, tc := range prefixTestCases() {
		t.Run(fmt.Sprintf("level_%d_rows_%d", tc.level, tc.rows), func(t *testing.T) {
			a, b := ring.NewPoly(params.N(), tc.level), ring.NewPoly(params.N(), tc.level)
			fillPrefixPoly(a, ringQ, 11)
			fillPrefixPoly(b, ringQ, 37)
			for _, sub := range []bool{false, true} {
				out := ring.NewPoly(params.N(), tc.level)
				fillPrefixSentinel(out, ringQ, 71)
				before := clonePrefixPoly(out)
				require.NoError(t, addSubPrefixRows(ringQ, tc.level, tc.rows, a, b, out, sub))
				for row := 0; row < tc.rows; row++ {
					want := make([]uint64, params.N())
					if sub {
						ringQ.SubRings[row].Sub(a.Coeffs[row], b.Coeffs[row], want)
					} else {
						ringQ.SubRings[row].Add(a.Coeffs[row], b.Coeffs[row], want)
					}
					require.Equal(t, want, out.Coeffs[row])
				}
				for row := tc.rows; row <= tc.level; row++ {
					require.Equal(t, before.Coeffs[row], out.Coeffs[row], "unrequested q%d", row)
				}
			}

			copied := ring.NewPoly(params.N(), tc.level)
			fillPrefixSentinel(copied, ringQ, 83)
			beforeCopy := clonePrefixPoly(copied)
			require.NoError(t, copyPrefixRows(ringQ, tc.level, tc.rows, a, copied))
			requirePrefixRowsEqual(t, a, copied, tc.rows)
			for row := tc.rows; row <= tc.level; row++ {
				require.Equal(t, beforeCopy.Coeffs[row], copied.Coeffs[row])
			}

			negated := ring.NewPoly(params.N(), tc.level)
			fillPrefixSentinel(negated, ringQ, 89)
			beforeNeg := clonePrefixPoly(negated)
			require.NoError(t, negatePrefixRows(ringQ, tc.level, tc.rows, a, negated))
			for row := 0; row < tc.rows; row++ {
				want := make([]uint64, params.N())
				ringQ.SubRings[row].Neg(a.Coeffs[row], want)
				require.Equal(t, want, negated.Coeffs[row])
			}
			for row := tc.rows; row <= tc.level; row++ {
				require.Equal(t, beforeNeg.Coeffs[row], negated.Coeffs[row])
			}

			in := NewCiphertext(params, 1, tc.level)
			outInteger := NewCiphertext(params, 1, tc.level)
			fillPrefixCiphertext(in, ringQ, 101)
			fillPrefixCiphertext(outInteger, ringQ, 131)
			beforeInteger := outInteger.CopyNew()
			integer := big.NewInt(-17)
			require.NoError(t, eval.mulIntegerRows(in, integer, outInteger, tc.rows))
			for component := range in.Value {
				for row := 0; row < tc.rows; row++ {
					q := new(big.Int).SetUint64(ringQ.SubRings[row].Modulus)
					for i, value := range in.Value[component].Coeffs[row] {
						want := new(big.Int).Mul(new(big.Int).SetUint64(value), integer)
						want.Mod(want, q)
						require.Equal(t, want.Uint64(), outInteger.Value[component].Coeffs[row][i])
					}
				}
				for row := tc.rows; row <= tc.level; row++ {
					require.Equal(t, beforeInteger.Value[component].Coeffs[row], outInteger.Value[component].Coeffs[row])
				}
			}
		})
	}
}

func TestPrefixCompletePointwiseMultiplyAndFusedAdd(t *testing.T) {
	params := testFastCKKSParameters(t)
	ringQ := params.RingQ()
	for _, tc := range prefixTestCases() {
		t.Run(fmt.Sprintf("level_%d_rows_%d", tc.level, tc.rows), func(t *testing.T) {
			a, b := ring.NewPoly(params.N(), tc.level), ring.NewPoly(params.N(), tc.level)
			out := ring.NewPoly(params.N(), tc.level)
			fillPrefixPoly(a, ringQ, 149)
			fillPrefixPoly(b, ringQ, 173)
			fillPrefixSentinel(out, ringQ, 197)
			before := clonePrefixPoly(out)
			require.NoError(t, pointMulPrefixRows(ringQ, tc.level, tc.rows, a, b, out, false, false))
			for row := 0; row < tc.rows; row++ {
				want := make([]uint64, params.N())
				ringQ.SubRings[row].MulCoeffsBarrett(a.Coeffs[row], b.Coeffs[row], want)
				require.Equal(t, want, out.Coeffs[row])
			}
			for row := tc.rows; row <= tc.level; row++ {
				require.Equal(t, before.Coeffs[row], out.Coeffs[row])
			}

			fused := clonePrefixPoly(before)
			require.NoError(t, pointMulPrefixRows(ringQ, tc.level, tc.rows, a, b, fused, false, true))
			for row := 0; row < tc.rows; row++ {
				want := clonePrefixPoly(before)
				ringQ.SubRings[row].MulCoeffsBarrettThenAdd(a.Coeffs[row], b.Coeffs[row], want.Coeffs[row])
				require.Equal(t, want.Coeffs[row], fused.Coeffs[row])
			}
			for row := tc.rows; row <= tc.level; row++ {
				require.Equal(t, before.Coeffs[row], fused.Coeffs[row])
			}
		})
	}
}

func TestPrefixCompletePointwiseMontgomeryRows(t *testing.T) {
	params := testFastCKKSParameters(t)
	ringQ := params.RingQ()
	for _, rows := range []int{2, 4} {
		t.Run(fmt.Sprintf("rows_%d", rows), func(t *testing.T) {
			const level = 3
			a, b := ring.NewPoly(params.N(), level), ring.NewPoly(params.N(), level)
			acc, out := ring.NewPoly(params.N(), level), ring.NewPoly(params.N(), level)
			fillPrefixPoly(a, ringQ, uint64(239+rows))
			fillPrefixPoly(b, ringQ, uint64(241+rows))
			fillPrefixPoly(acc, ringQ, uint64(251+rows))
			fillPrefixSentinel(out, ringQ, uint64(257+rows))
			for row := 0; row <= level; row++ {
				ringQ.SubRings[row].MForm(a.Coeffs[row], a.Coeffs[row])
				ringQ.SubRings[row].MForm(b.Coeffs[row], b.Coeffs[row])
				ringQ.SubRings[row].MForm(acc.Coeffs[row], acc.Coeffs[row])
			}
			beforeAcc, beforeOut := clonePrefixPoly(acc), clonePrefixPoly(out)
			require.NoError(t, pointMulPrefixRows(ringQ, level, rows, a, b, out, true, false))
			require.NoError(t, pointMulPrefixRows(ringQ, level, rows, a, b, acc, true, true))
			for row := 0; row < rows; row++ {
				wantProduct := make([]uint64, params.N())
				ringQ.SubRings[row].MulCoeffsMontgomery(a.Coeffs[row], b.Coeffs[row], wantProduct)
				wantAccum := append([]uint64(nil), beforeAcc.Coeffs[row]...)
				ringQ.SubRings[row].MulCoeffsMontgomeryThenAdd(a.Coeffs[row], b.Coeffs[row], wantAccum)
				require.Equal(t, wantProduct, out.Coeffs[row])
				require.Equal(t, wantAccum, acc.Coeffs[row])
			}
			for row := rows; row <= level; row++ {
				require.Equal(t, beforeOut.Coeffs[row], out.Coeffs[row])
				require.Equal(t, beforeAcc.Coeffs[row], acc.Coeffs[row])
			}
		})
	}
}

func TestPrefixCompleteDegreeOneMulMontgomeryRows(t *testing.T) {
	params := testFastCKKSParameters(t)
	ringQ := params.RingQ()
	eval := NewEvaluator(params)
	const level = 3
	a, b := NewCiphertext(params, 1, level), NewCiphertext(params, 1, level)
	out, relinOut := NewCiphertext(params, 2, level), NewCiphertext(params, 1, level)
	fillPrefixCiphertext(a, ringQ, 263)
	fillPrefixCiphertext(b, ringQ, 269)
	for _, ct := range []*rlwe.Ciphertext{a, b} {
		for component := range ct.Value {
			for row := 0; row <= level; row++ {
				ringQ.SubRings[row].MForm(ct.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
			}
		}
		ct.IsNTT, ct.IsMontgomery = true, true
	}
	out.IsNTT, out.IsMontgomery = true, true
	relinOut.IsNTT, relinOut.IsMontgomery = true, true
	require.NoError(t, eval.mulElementRows(a, b.El(), out, false, MaxQPrefixWidth))
	require.NoError(t, eval.mulElementRows(a, b.El(), relinOut, true, MaxQPrefixWidth))
	for row := 0; row < MaxQPrefixWidth; row++ {
		subring := ringQ.SubRings[row]
		want0, want1a, want1b, want2 := make([]uint64, params.N()), make([]uint64, params.N()), make([]uint64, params.N()), make([]uint64, params.N())
		subring.MulCoeffsMontgomery(a.Value[0].Coeffs[row], b.Value[0].Coeffs[row], want0)
		subring.MulCoeffsMontgomery(a.Value[0].Coeffs[row], b.Value[1].Coeffs[row], want1a)
		subring.MulCoeffsMontgomery(a.Value[1].Coeffs[row], b.Value[0].Coeffs[row], want1b)
		subring.Add(want1a, want1b, want1a)
		subring.MulCoeffsMontgomery(a.Value[1].Coeffs[row], b.Value[1].Coeffs[row], want2)
		require.Equal(t, want0, out.Value[0].Coeffs[row])
		require.Equal(t, want1a, out.Value[1].Coeffs[row])
		require.Equal(t, want2, out.Value[2].Coeffs[row])
		require.Equal(t, want0, relinOut.Value[0].Coeffs[row])
		require.Equal(t, want1a, relinOut.Value[1].Coeffs[row])
	}
}

func TestPrefixCompletePartialNTTINTTAndAutomorphism(t *testing.T) {
	params := testFastCKKSParameters(t)
	ringQ := params.RingQ()
	eval := NewEvaluator(params)
	galEl := params.GaloisElement(1)
	for _, tc := range prefixTestCases() {
		t.Run(fmt.Sprintf("level_%d_rows_%d", tc.level, tc.rows), func(t *testing.T) {
			coeff := ring.NewPoly(params.N(), tc.level)
			fillPrefixPoly(coeff, ringQ, 211)
			ntt := ring.NewPoly(params.N(), tc.level)
			fillPrefixSentinel(ntt, ringQ, 223)
			beforeNTT := clonePrefixPoly(ntt)
			require.NoError(t, FastPartialNTTRows(ringQ, coeff, ntt, tc.rows))
			wantNTT := clonePrefixPoly(coeff)
			ringQ.AtLevel(tc.level).NTT(coeff, wantNTT)
			requirePrefixRowsEqual(t, wantNTT, ntt, tc.rows)
			for row := tc.rows; row <= tc.level; row++ {
				require.Equal(t, beforeNTT.Coeffs[row], ntt.Coeffs[row])
			}

			coeffOut := ring.NewPoly(params.N(), tc.level)
			fillPrefixSentinel(coeffOut, ringQ, 227)
			beforeINTT := clonePrefixPoly(coeffOut)
			require.NoError(t, FastPartialINTTRows(ringQ, ntt, coeffOut, tc.rows))
			wantCoeff := clonePrefixPoly(ntt)
			ringQ.AtLevel(tc.level).INTT(ntt, wantCoeff)
			requirePrefixRowsEqual(t, wantCoeff, coeffOut, tc.rows)
			for row := tc.rows; row <= tc.level; row++ {
				require.Equal(t, beforeINTT.Coeffs[row], coeffOut.Coeffs[row])
			}

			auto := ring.NewPoly(params.N(), tc.level)
			fillPrefixSentinel(auto, ringQ, 229)
			beforeAuto := clonePrefixPoly(auto)
			require.NoError(t, FastAutomorphismRows(ringQ, ntt, auto, galEl, true, tc.rows))
			wantAuto := clonePrefixPoly(ntt)
			ringQ.AtLevel(tc.level).AutomorphismNTT(ntt, galEl, wantAuto)
			requirePrefixRowsEqual(t, wantAuto, auto, tc.rows)
			for row := tc.rows; row <= tc.level; row++ {
				require.Equal(t, beforeAuto.Coeffs[row], auto.Coeffs[row])
			}

			scratchAuto := ring.NewPoly(params.N(), tc.level)
			fillPrefixSentinel(scratchAuto, ringQ, 233)
			require.NoError(t, eval.fastAutomorphismRows(ringQ, ntt, scratchAuto, galEl, true, tc.rows))
			requirePrefixRowsEqual(t, wantAuto, scratchAuto, tc.rows)
		})
	}
}

func TestPrefixCompleteScalarNTTEncodingUsesActualModuli(t *testing.T) {
	params := testFastCKKSParameters(t)
	eval := NewEvaluator(params)
	scalar := bignum.ToComplex(complex(1.25, -0.5), params.EncodingPrecision())
	scale := new(big.Float).SetInt64(8)
	for _, tc := range prefixTestCases() {
		for _, montgomery := range []bool{false, true} {
			t.Run(fmt.Sprintf("level_%d_rows_%d_montgomery_%t", tc.level, tc.rows, montgomery), func(t *testing.T) {
				got, err := eval.scalarNTTRows(scalar, scale, montgomery, tc.level, tc.rows)
				require.NoError(t, err)
				real := big.NewInt(10)
				imag := big.NewInt(-4)
				for row := 0; row < tc.rows; row++ {
					subring := params.RingQ().SubRings[row]
					modulus := new(big.Int).SetUint64(subring.Modulus)
					realMod := new(big.Int).Mod(new(big.Int).Set(real), modulus).Uint64()
					imagMod := new(big.Int).Mod(new(big.Int).Set(imag), modulus).Uint64()
					imagMont := ring.MRed(imagMod, subring.RootsForward[1], subring.Modulus, subring.MRedConstant)
					want0 := ring.CRed(realMod+imagMont, subring.Modulus)
					want1 := ring.CRed(realMod+subring.Modulus-imagMont, subring.Modulus)
					if montgomery {
						want0 = ring.MForm(want0, subring.Modulus, subring.BRedConstant)
						want1 = ring.MForm(want1, subring.Modulus, subring.BRedConstant)
					}
					require.Equal(t, want0, got[row][0], "q%d positive half", row)
					require.Equal(t, want1, got[row][1], "q%d negative half", row)
				}
				for row := tc.rows; row < MaxQPrefixWidth; row++ {
					require.Zero(t, got[row])
				}
			})
		}
	}
}

func TestPrefixCompleteScalarMultiplyRows(t *testing.T) {
	params := testFastCKKSParameters(t)
	ringQ := params.RingQ()
	eval := NewEvaluator(params)
	scalar := bignum.ToComplex(uint64(3), params.EncodingPrecision())
	scale := rlwe.NewScale(1)
	for _, tc := range prefixTestCases() {
		t.Run(fmt.Sprintf("level_%d_rows_%d", tc.level, tc.rows), func(t *testing.T) {
			in, out := NewCiphertext(params, 1, tc.level), NewCiphertext(params, 1, tc.level)
			fillPrefixCiphertext(in, ringQ, uint64(487+tc.rows))
			fillPrefixCiphertext(out, ringQ, uint64(491+tc.rows))
			in.IsNTT, out.IsNTT = true, true
			before := out.CopyNew()
			require.NoError(t, eval.mulScalarAtScaleRows(in, scalar, scale, out, tc.rows))
			for component := range in.Value {
				for row := 0; row < tc.rows; row++ {
					q := new(big.Int).SetUint64(ringQ.SubRings[row].Modulus)
					for i, value := range in.Value[component].Coeffs[row] {
						wantBig := new(big.Int).Mul(new(big.Int).SetUint64(value), big.NewInt(3))
						want := wantBig.Mod(wantBig, q).Uint64()
						require.Equal(t, want, out.Value[component].Coeffs[row][i])
					}
				}
				for row := tc.rows; row <= tc.level; row++ {
					require.Equal(t, before.Value[component].Coeffs[row], out.Value[component].Coeffs[row])
				}
			}
		})
	}
}

func TestPrefixCompleteTruncateAndDegreeOneMulFormulas(t *testing.T) {
	params := testFastCKKSParameters(t)
	ringQ := params.RingQ()
	eval := NewEvaluator(params)

	for _, tc := range prefixTestCases() {
		t.Run(fmt.Sprintf("truncate_level_%d_rows_%d", tc.level, tc.rows), func(t *testing.T) {
			in := NewCiphertext(params, 2, tc.level)
			out := NewCiphertext(params, 1, tc.level)
			fillPrefixCiphertext(in, ringQ, 251)
			fillPrefixCiphertext(out, ringQ, 271)
			before := out.CopyNew()
			require.NoError(t, fastTruncateDegree2To1Rows(ringQ, in, out, tc.rows))
			require.Equal(t, 1, out.Degree())
			for component := 0; component < 2; component++ {
				for row := 0; row < tc.rows; row++ {
					require.Equal(t, in.Value[component].Coeffs[row], out.Value[component].Coeffs[row])
				}
				for row := tc.rows; row <= tc.level; row++ {
					require.Equal(t, before.Value[component].Coeffs[row], out.Value[component].Coeffs[row])
				}
			}
		})
	}

	for _, tc := range prefixTestCases() {
		t.Run(fmt.Sprintf("mul_level_%d_rows_%d", tc.level, tc.rows), func(t *testing.T) {
			a, b := NewCiphertext(params, 1, tc.level), NewCiphertext(params, 1, tc.level)
			out, relinOut := NewCiphertext(params, 2, tc.level), NewCiphertext(params, 1, tc.level)
			accum := NewCiphertext(params, 2, tc.level)
			fillPrefixCiphertext(a, ringQ, uint64(281+tc.rows))
			fillPrefixCiphertext(b, ringQ, uint64(307+tc.rows))
			fillPrefixCiphertext(out, ringQ, uint64(331+tc.rows))
			fillPrefixCiphertext(relinOut, ringQ, uint64(353+tc.rows))
			fillPrefixCiphertext(accum, ringQ, uint64(359+tc.rows))
			a.IsNTT, b.IsNTT, out.IsNTT, relinOut.IsNTT, accum.IsNTT = true, true, true, true, true
			beforeOut, beforeRelin, beforeAccum := out.CopyNew(), relinOut.CopyNew(), accum.CopyNew()
			require.NoError(t, eval.mulElementRows(a, b.El(), out, false, tc.rows))
			require.NoError(t, eval.mulElementRows(a, b.El(), relinOut, true, tc.rows))
			require.NoError(t, eval.mulElementThenAdd(a, b.El(), accum, tc.rows))
			require.Equal(t, 2, out.Degree())
			require.Equal(t, 1, relinOut.Degree())
			for row := 0; row < tc.rows; row++ {
				q := new(big.Int).SetUint64(ringQ.SubRings[row].Modulus)
				for i := 0; i < params.N(); i++ {
					qMod := ringQ.SubRings[row].Modulus
					want0 := prefixMulMod(a.Value[0].Coeffs[row][i], b.Value[0].Coeffs[row][i], qMod)
					want1a := prefixMulMod(a.Value[0].Coeffs[row][i], b.Value[1].Coeffs[row][i], qMod)
					want1b := prefixMulMod(a.Value[1].Coeffs[row][i], b.Value[0].Coeffs[row][i], qMod)
					want1Big := new(big.Int).Add(new(big.Int).SetUint64(want1a), new(big.Int).SetUint64(want1b))
					want1 := want1Big.Mod(want1Big, q).Uint64()
					want2 := prefixMulMod(a.Value[1].Coeffs[row][i], b.Value[1].Coeffs[row][i], qMod)
					require.Equal(t, want0, out.Value[0].Coeffs[row][i])
					require.Equal(t, want1, out.Value[1].Coeffs[row][i])
					require.Equal(t, want2, out.Value[2].Coeffs[row][i])
					require.Equal(t, want0, relinOut.Value[0].Coeffs[row][i])
					require.Equal(t, want1, relinOut.Value[1].Coeffs[row][i])
					require.Equal(t, prefixAddMod(beforeAccum.Value[0].Coeffs[row][i], want0, qMod), accum.Value[0].Coeffs[row][i])
					require.Equal(t, prefixAddMod(beforeAccum.Value[1].Coeffs[row][i], want1, qMod), accum.Value[1].Coeffs[row][i])
					require.Equal(t, prefixAddMod(beforeAccum.Value[2].Coeffs[row][i], want2, qMod), accum.Value[2].Coeffs[row][i])
				}
			}
			for component := range out.Value {
				for row := tc.rows; row <= tc.level; row++ {
					require.Equal(t, beforeOut.Value[component].Coeffs[row], out.Value[component].Coeffs[row])
				}
			}
			for component := range relinOut.Value {
				for row := tc.rows; row <= tc.level; row++ {
					require.Equal(t, beforeRelin.Value[component].Coeffs[row], relinOut.Value[component].Coeffs[row])
				}
			}
			for component := range accum.Value {
				for row := tc.rows; row <= tc.level; row++ {
					require.Equal(t, beforeAccum.Value[component].Coeffs[row], accum.Value[component].Coeffs[row])
				}
			}
		})
	}
}

func TestPrefixRowValidationRejectsInvalidWidthBeforeWriting(t *testing.T) {
	params := testFastCKKSParameters(t)
	ringQ := params.RingQ()
	a, b, out := ring.NewPoly(params.N(), 3), ring.NewPoly(params.N(), 3), ring.NewPoly(params.N(), 3)
	fillPrefixPoly(a, ringQ, 367)
	fillPrefixPoly(b, ringQ, 373)
	fillPrefixSentinel(out, ringQ, 379)
	before := clonePrefixPoly(out)
	require.Error(t, pointMulPrefixRows(ringQ, 3, 5, a, b, out, false, false))
	require.Error(t, pointMulPrefixRows(ringQ, 3, 0, a, b, out, false, false))
	require.Equal(t, before, out)

	missing := clonePrefixPoly(a)
	missing.Coeffs[2] = nil
	require.Error(t, FastPartialNTTRows(ringQ, missing, out, 3))
	require.Equal(t, before, out)
}

func TestFastLegacyWrappersDoNotPromotePoisonedQ3(t *testing.T) {
	q0Generator := ring.NewNTTFriendlyPrimesGenerator(55, 1<<5)
	q0, err := q0Generator.NextUpstreamPrime()
	require.NoError(t, err)
	q123Generator := ring.NewNTTFriendlyPrimesGenerator(39, 1<<5)
	q123, err := q123Generator.NextDownstreamPrimes(3)
	require.NoError(t, err)
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            4,
		Q:               append([]uint64{q0}, q123...),
		LogDefaultScale: 30,
	})
	require.NoError(t, err)
	require.True(t, q012Enabled(params), "Q=%v", params.Q())
	require.Equal(t, 3, maintainedLimbCount(&params, 3))
	eval := NewEvaluator(params)
	ringQ := params.RingQ()

	a, b := NewCiphertext(params, 1, 3), NewCiphertext(params, 1, 3)
	fillPrefixCiphertext(a, ringQ, 401)
	fillPrefixCiphertext(b, ringQ, 431)
	for i := range a.Value[0].Coeffs[3] {
		a.Value[0].Coeffs[3][i] = uint64(0xabc000 + i)
		b.Value[0].Coeffs[3][i] = uint64(0xdef000 + i)
	}
	a.IsNTT, b.IsNTT = true, true
	for _, ct := range []*rlwe.Ciphertext{a, b} {
		ringQ.AtLevel(3).NTT(ct.Value[0], ct.Value[0])
		ringQ.AtLevel(3).NTT(ct.Value[1], ct.Value[1])
	}
	for i := range a.Value[0].Coeffs[3] {
		a.Value[0].Coeffs[3][i] = uint64(0xabc000 + i)
		b.Value[0].Coeffs[3][i] = uint64(0xdef000 + i)
	}

	addOut := NewCiphertext(params, 1, 3)
	addOut.IsNTT = true
	addQ3 := make([][]uint64, len(addOut.Value))
	for component := range addOut.Value {
		for i := range addOut.Value[component].Coeffs[3] {
			addOut.Value[component].Coeffs[3][i] = uint64(0x123000 + 100*component + i)
		}
		addQ3[component] = append([]uint64(nil), addOut.Value[component].Coeffs[3]...)
	}
	require.NoError(t, eval.Add(a, b, addOut))
	for row := 0; row < 3; row++ {
		for component := range addOut.Value {
			want := make([]uint64, params.N())
			ringQ.SubRings[row].Add(a.Value[component].Coeffs[row], b.Value[component].Coeffs[row], want)
			require.Equal(t, want, addOut.Value[component].Coeffs[row])
		}
	}
	for component := range addOut.Value {
		require.Equal(t, addQ3[component], addOut.Value[component].Coeffs[3])
	}

	mulOut := NewCiphertext(params, 2, 3)
	mulOut.IsNTT = true
	for component := range mulOut.Value {
		for i := range mulOut.Value[component].Coeffs[3] {
			mulOut.Value[component].Coeffs[3][i] = uint64(0x234000 + 100*component + i)
		}
	}
	mulQ3 := make([][]uint64, len(mulOut.Value))
	for component := range mulOut.Value {
		mulQ3[component] = append([]uint64(nil), mulOut.Value[component].Coeffs[3]...)
	}
	require.NoError(t, eval.Mul(a, b, mulOut))
	for component := range mulOut.Value {
		require.Equal(t, mulQ3[component], mulOut.Value[component].Coeffs[3])
	}
	for row := 0; row < 3; row++ {
		qMod := ringQ.SubRings[row].Modulus
		for i := 0; i < params.N(); i++ {
			want0 := prefixMulMod(a.Value[0].Coeffs[row][i], b.Value[0].Coeffs[row][i], qMod)
			want1a := prefixMulMod(a.Value[0].Coeffs[row][i], b.Value[1].Coeffs[row][i], qMod)
			want1b := prefixMulMod(a.Value[1].Coeffs[row][i], b.Value[0].Coeffs[row][i], qMod)
			want1Big := new(big.Int).Add(new(big.Int).SetUint64(want1a), new(big.Int).SetUint64(want1b))
			want1 := want1Big.Mod(want1Big, new(big.Int).SetUint64(qMod)).Uint64()
			want2 := prefixMulMod(a.Value[1].Coeffs[row][i], b.Value[1].Coeffs[row][i], qMod)
			require.Equal(t, want0, mulOut.Value[0].Coeffs[row][i])
			require.Equal(t, want1, mulOut.Value[1].Coeffs[row][i])
			require.Equal(t, want2, mulOut.Value[2].Coeffs[row][i])
		}
	}

	scalarOut := NewCiphertext(params, 1, 3)
	scalarOut.IsNTT = true
	scalarQ3 := make([][]uint64, len(scalarOut.Value))
	for component := range scalarOut.Value {
		for i := range scalarOut.Value[component].Coeffs[3] {
			scalarOut.Value[component].Coeffs[3][i] = uint64(0x345000 + 100*component + i)
		}
		scalarQ3[component] = append([]uint64(nil), scalarOut.Value[component].Coeffs[3]...)
	}
	require.NoError(t, eval.Mul(a, uint64(3), scalarOut))
	for component := range scalarOut.Value {
		require.Equal(t, scalarQ3[component], scalarOut.Value[component].Coeffs[3])
		for row := 0; row < 3; row++ {
			q := new(big.Int).SetUint64(ringQ.SubRings[row].Modulus)
			for i, value := range a.Value[component].Coeffs[row] {
				want := new(big.Int).Mul(new(big.Int).SetUint64(value), big.NewInt(3))
				want.Mod(want, q)
				require.Equal(t, want.Uint64(), scalarOut.Value[component].Coeffs[row][i])
			}
		}
	}

	partialIn, partialOut := ring.NewPoly(params.N(), 3), ring.NewPoly(params.N(), 3)
	fillPrefixPoly(partialIn, ringQ, 443)
	fillPrefixSentinel(partialOut, ringQ, 449)
	partialQ3 := append([]uint64(nil), partialOut.Coeffs[3]...)
	require.NoError(t, FastPartialNTT(ringQ, partialIn, partialOut))
	wantPartial := clonePrefixPoly(partialIn)
	ringQ.NTT(partialIn, wantPartial)
	requirePrefixRowsEqual(t, wantPartial, partialOut, 3)
	require.Equal(t, partialQ3, partialOut.Coeffs[3])

	autoOut := ring.NewPoly(params.N(), 3)
	fillPrefixSentinel(autoOut, ringQ, 457)
	autoQ3 := append([]uint64(nil), autoOut.Coeffs[3]...)
	require.NoError(t, FastAutomorphism(ringQ, partialIn, autoOut, params.GaloisElement(1), false))
	wantAuto := clonePrefixPoly(partialIn)
	ringQ.Automorphism(partialIn, params.GaloisElement(1), wantAuto)
	requirePrefixRowsEqual(t, wantAuto, autoOut, 3)
	require.Equal(t, autoQ3, autoOut.Coeffs[3])

	degreeTwo, truncOut := NewCiphertext(params, 2, 3), NewCiphertext(params, 1, 3)
	fillPrefixCiphertext(degreeTwo, ringQ, 461)
	fillPrefixCiphertext(truncOut, ringQ, 467)
	degreeTwo.IsNTT, truncOut.IsNTT = true, true
	truncQ3 := make([][]uint64, len(truncOut.Value))
	for component := range truncOut.Value {
		truncQ3[component] = append([]uint64(nil), truncOut.Value[component].Coeffs[3]...)
	}
	require.NoError(t, eval.Relinearize(degreeTwo, truncOut))
	for component := range truncOut.Value {
		require.Equal(t, truncQ3[component], truncOut.Value[component].Coeffs[3])
		for row := 0; row < 3; row++ {
			require.Equal(t, degreeTwo.Value[component].Coeffs[row], truncOut.Value[component].Coeffs[row])
		}
	}
}
