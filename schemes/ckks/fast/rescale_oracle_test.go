package fast

import (
	"fmt"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestFastRescaleMatchesBigIntCenteredCRTOracleAtEveryPrefixWidth(t *testing.T) {
	params := rescaleHighLevelParameters(t)
	fastEval := NewEvaluator(params)
	standardEval := ckks.NewEvaluator(params, nil)
	for level := 1; level <= 4; level++ {
		t.Run(fmt.Sprintf("level_%d", level), func(t *testing.T) {
			sourceRows, err := QPrefixWidth(level)
			require.NoError(t, err)
			values := rescaleOracleBoundaryValues(t, params, level, sourceRows)
			fastIn, standardIn := makeRescaleOracleInputs(params, level, sourceRows, values)
			fastOut := NewCiphertext(params, fastIn.Degree(), level-1)
			standardOut := ckks.NewCiphertext(params, standardIn.Degree(), level-1)
			require.NoError(t, fastEval.RescaleQPrefixRows(fastIn, sourceRows, fastOut))
			require.NoError(t, standardEval.Rescale(standardIn, standardOut))
			require.Equal(t, level-1, fastOut.Level())
			require.Equal(t, standardOut.Level(), fastOut.Level())
			require.Equal(t, standardOut.Scale, fastOut.Scale)
			want := bigIntRescaleOracle(t, params, fastIn, 1, sourceRows)
			for component := range fastOut.Value {
				for row := range want[component] {
					require.Equal(t, want[component][row], fastOut.Value[component].Coeffs[row], "oracle component=%d q%d", component, row)
					require.Equal(t, standardOut.Value[component].Coeffs[row], fastOut.Value[component].Coeffs[row], "Standard component=%d q%d", component, row)
				}
			}
		})
	}
}

func TestFastRescaleToMatchesSequentialBigIntOracleAcrossContractions(t *testing.T) {
	params := rescaleHighLevelParameters(t)
	fastEval := NewEvaluator(params)
	standardEval := ckks.NewEvaluator(params, nil)
	for _, tc := range []struct{ level, steps int }{{4, 2}, {3, 2}, {2, 2}} {
		for _, authority := range []struct {
			name     string
			rows     int
			explicit bool
		}{{"default_qprefix", qPrefixWidthOrPanic(tc.level), false}, {"explicit_qprefix", qPrefixWidthOrPanic(tc.level), true}} {
			t.Run(fmt.Sprintf("level_%d_steps_%d/%s", tc.level, tc.steps, authority.name), func(t *testing.T) {
				values := rescaleOracleBoundaryValues(t, params, tc.level, authority.rows)
				fastIn, standardIn := makeRescaleOracleInputs(params, tc.level, authority.rows, values)
				scaleValue := new(big.Int).Lsh(big.NewInt(1), 30)
				for step := 0; step < tc.steps; step++ {
					scaleValue.Mul(scaleValue, new(big.Int).SetUint64(params.Q()[tc.level-step]))
				}
				fastIn.Scale, standardIn.Scale = rlwe.NewScale(scaleValue), rlwe.NewScale(scaleValue)
				minScale := rlwe.NewScale(1 << 30)
				targetLevel := tc.level - tc.steps
				fastOut := NewCiphertext(params, fastIn.Degree(), tc.level)
				standardOut := ckks.NewCiphertext(params, standardIn.Degree(), tc.level)
				var err error
				if authority.explicit {
					err = fastEval.RescaleToQPrefixRows(fastIn, minScale, authority.rows, fastOut)
				} else {
					err = fastEval.RescaleTo(fastIn, minScale, fastOut)
				}
				require.NoError(t, err)
				require.NoError(t, standardEval.RescaleTo(standardIn, minScale, standardOut))
				require.Equal(t, targetLevel, fastOut.Level())
				require.Equal(t, standardOut.Level(), fastOut.Level())
				require.Equal(t, standardOut.Scale, fastOut.Scale)
				want := bigIntRescaleOracle(t, params, fastIn, tc.steps, authority.rows)
				for component := range fastOut.Value {
					for row := range want[component] {
						require.Equal(t, want[component][row], fastOut.Value[component].Coeffs[row], "oracle component=%d q%d", component, row)
						require.Equal(t, standardOut.Value[component].Coeffs[row], fastOut.Value[component].Coeffs[row], "Standard component=%d q%d", component, row)
					}
				}
			})
		}
	}
}

func rescaleOracleBoundaryValues(t testing.TB, params ckks.Parameters, level, sourceRows int) []*big.Int {
	t.Helper()
	sourceProduct := big.NewInt(1)
	for row := 0; row < sourceRows; row++ {
		sourceProduct.Mul(sourceProduct, new(big.Int).SetUint64(params.Q()[row]))
	}
	sourceHalf := new(big.Int).Rsh(new(big.Int).Set(sourceProduct), 1)
	divisor := params.Q()[level]
	divisorHalf := new(big.Int).SetUint64(divisor / 2)
	values := []*big.Int{
		big.NewInt(0), big.NewInt(1), big.NewInt(-1),
		new(big.Int).Sub(new(big.Int).Set(divisorHalf), big.NewInt(1)),
		new(big.Int).Set(divisorHalf),
		new(big.Int).Add(new(big.Int).Set(divisorHalf), big.NewInt(1)),
		new(big.Int).Neg(new(big.Int).Set(divisorHalf)),
		new(big.Int).Neg(new(big.Int).Add(new(big.Int).Set(divisorHalf), big.NewInt(1))),
		new(big.Int).Sub(new(big.Int).Set(sourceHalf), big.NewInt(1)),
		new(big.Int).Set(sourceHalf),
		new(big.Int).Add(new(big.Int).Set(sourceHalf), big.NewInt(1)),
		new(big.Int).Neg(new(big.Int).Set(sourceHalf)),
	}
	if level < 4 {
		targetRows := min(sourceRows, qPrefixWidthOrPanic(level-1))
		targetProduct := big.NewInt(1)
		for row := 0; row < targetRows; row++ {
			targetProduct.Mul(targetProduct, new(big.Int).SetUint64(params.Q()[row]))
		}
		targetHalf := new(big.Int).Rsh(new(big.Int).Set(targetProduct), 1)
		nearCapacity := new(big.Int).Mul(targetHalf, new(big.Int).SetUint64(divisor))
		nearCapacity.Add(nearCapacity, divisorHalf)
		for _, candidate := range []*big.Int{nearCapacity, new(big.Int).Sub(new(big.Int).Set(nearCapacity), big.NewInt(1))} {
			values = append(values, candidate, new(big.Int).Neg(new(big.Int).Set(candidate)))
		}
	}
	if len(values) > params.N() {
		values = values[:params.N()]
	}
	for len(values) < params.N() {
		values = append(values, new(big.Int))
	}
	return values
}

func makeRescaleOracleInputs(params ckks.Parameters, level, sourceRows int, values []*big.Int) (fast, standard *rlwe.Ciphertext) {
	fast = NewCiphertext(params, 1, level)
	standard = ckks.NewCiphertext(params, 1, level)
	fast.IsNTT, standard.IsNTT = true, true
	fast.Scale, standard.Scale = rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 180)), rlwe.NewScale(new(big.Int).Lsh(big.NewInt(1), 180))
	for component := range fast.Value {
		for row := 0; row <= level; row++ {
			q := params.Q()[row]
			for coefficient := 0; coefficient < params.N(); coefficient++ {
				signed := new(big.Int).Set(values[coefficient])
				if component == 1 {
					signed.Neg(signed)
				}
				canonical := centeredValueForPrefix(signed, params.Q()[:sourceRows])
				residue := new(big.Int).Mod(canonical, new(big.Int).SetUint64(q)).Uint64()
				standard.Value[component].Coeffs[row][coefficient] = residue
				if row < sourceRows {
					fast.Value[component].Coeffs[row][coefficient] = residue
				}
			}
			params.RingQ().SubRings[row].NTT(standard.Value[component].Coeffs[row], standard.Value[component].Coeffs[row])
			if row < sourceRows {
				params.RingQ().SubRings[row].NTT(fast.Value[component].Coeffs[row], fast.Value[component].Coeffs[row])
			}
		}
	}
	return fast, standard
}

func centeredValueForPrefix(value *big.Int, q []uint64) *big.Int {
	product := big.NewInt(1)
	for _, modulus := range q {
		product.Mul(product, new(big.Int).SetUint64(modulus))
	}
	canonical := new(big.Int).Mod(new(big.Int).Set(value), product)
	if canonical.Cmp(new(big.Int).Rsh(new(big.Int).Set(product), 1)) > 0 {
		canonical.Sub(canonical, product)
	}
	return canonical
}

func bigIntRescaleOracle(t *testing.T, params ckks.Parameters, input *rlwe.Ciphertext, steps, sourceRows int) [][][]uint64 {
	t.Helper()
	level := input.Level()
	targetLevel := level - steps
	targetRows := min(sourceRows, qPrefixWidthOrPanic(targetLevel))
	ringQ := params.RingQ()
	want := make([][][]uint64, len(input.Value))
	for component := range input.Value {
		coeffRows := make([][]uint64, sourceRows)
		for row := 0; row < sourceRows; row++ {
			coeffRows[row] = make([]uint64, params.N())
			ringQ.SubRings[row].INTT(input.Value[component].Coeffs[row], coeffRows[row])
		}
		want[component] = make([][]uint64, targetRows)
		for row := 0; row < targetRows; row++ {
			want[component][row] = make([]uint64, params.N())
		}
		for coefficient := 0; coefficient < params.N(); coefficient++ {
			residues := make([]uint64, sourceRows)
			for row := 0; row < sourceRows; row++ {
				residues[row] = coeffRows[row][coefficient]
			}
			value := centeredBigIntCRT(residues, params.Q()[:sourceRows])
			for step := 0; step < steps; step++ {
				value = roundSignedBigInt(value, params.Q()[level-step])
				stepTarget := level - step - 1
				stepTargetRows := min(sourceRows, qPrefixWidthOrPanic(stepTarget))
				capacity := big.NewInt(1)
				for row := 0; row < stepTargetRows; row++ {
					capacity.Mul(capacity, new(big.Int).SetUint64(params.Q()[row]))
				}
				twiceMagnitude := new(big.Int).Lsh(new(big.Int).Abs(new(big.Int).Set(value)), 1)
				require.Less(t, twiceMagnitude.Cmp(capacity), 0, "exact result must fit Level %d capacity", stepTarget)
			}
			for row := 0; row < targetRows; row++ {
				want[component][row][coefficient] = new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(params.Q()[row])).Uint64()
			}
		}
		for row := 0; row < targetRows; row++ {
			ringQ.SubRings[row].NTT(want[component][row], want[component][row])
		}
	}
	return want
}

func centeredBigIntCRT(residues, q []uint64) *big.Int {
	x, product := new(big.Int), big.NewInt(1)
	for row, modulus := range q {
		qi := new(big.Int).SetUint64(modulus)
		xMod := new(big.Int).Mod(new(big.Int).Set(x), qi)
		delta := new(big.Int).Sub(new(big.Int).SetUint64(residues[row]), xMod)
		delta.Mod(delta, qi)
		productMod := new(big.Int).Mod(new(big.Int).Set(product), qi)
		inverse := new(big.Int).ModInverse(productMod, qi)
		if inverse == nil {
			panic("test Q-prefix moduli are not pairwise coprime")
		}
		digit := delta.Mul(delta, inverse)
		digit.Mod(digit, qi)
		x.Add(x, new(big.Int).Mul(product, digit))
		product.Mul(product, qi)
	}
	if x.Cmp(new(big.Int).Rsh(new(big.Int).Set(product), 1)) > 0 {
		x.Sub(x, product)
	}
	return x
}

func roundSignedBigInt(value *big.Int, divisor uint64) *big.Int {
	negative := value.Sign() < 0
	magnitude := new(big.Int).Abs(new(big.Int).Set(value))
	quotient, remainder := new(big.Int).QuoRem(magnitude, new(big.Int).SetUint64(divisor), new(big.Int))
	if new(big.Int).Lsh(remainder, 1).Cmp(new(big.Int).SetUint64(divisor)) > 0 {
		quotient.Add(quotient, big.NewInt(1))
	}
	if negative {
		quotient.Neg(quotient)
	}
	return quotient
}
