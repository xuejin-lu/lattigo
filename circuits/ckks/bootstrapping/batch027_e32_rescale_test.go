package bootstrapping

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	fastckks "github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

const (
	batch027FixtureName      = "batch027-e32-rescale-synthetic-v1"
	batch027FixtureSeed      = uint64(0x027e32c0ffee1234)
	batch027ConfigSHA256     = "919a2d9409b8ddeb720458d3112855f87d7a69e0cc39aad825b0ddf794769c98"
	batch027ExpectedQPSHA256 = "1f045e603a856968779d62e045a037274bba08cbfce8b1dd3dec2828f1f6a46b"
	batch027FixtureSHA256    = "1de1a696fe7b596105cf236e7799b015049ae78e63e09e2ad6314be730dc3cb8"
	batch027SourceLevel      = 9
	batch027SourcePrefixRows = 4
	batch027TargetLevel      = 8
	batch027TargetRows       = 4
)

type batch027SplitMix64 struct{ state uint64 }

func (rng *batch027SplitMix64) next() uint64 {
	// SplitMix64 v1; uint64 addition and multiplication wrap modulo 2^64.
	rng.state += 0x9e3779b97f4a7c15
	z := rng.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}

func batch027E32Parameters(tb testing.TB) (Parameters, ckks.Parameters) {
	tb.Helper()
	residual, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39}, LogDefaultScale: 45, Xs: ring.Ternary{H: 192},
	})
	require.NoError(tb, err)
	logN, logSlots, ephemeral, evalModScale, degree, doubleAngle := 13, 12, 32, 60, 30, 3
	k, logMessageRatio, mod1InvDegree := 16, 10, 0
	params, err := NewParametersFromLiteral(residual, ParametersLiteral{
		LogN: &logN, LogSlots: &logSlots, LogP: []int{61, 61, 61, 61, 61},
		SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{39}, {39}, {39}},
		CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{56}, {56}, {56}, {56}},
		EvalModLogScale: &evalModScale, EphemeralSecretWeight: &ephemeral, Mod1Type: mod1.CosDiscrete,
		LogMessageRatio: &logMessageRatio, K: &k, Mod1Degree: &degree, DoubleAngle: &doubleAngle,
		Mod1InvDegree: &mod1InvDegree,
	})
	require.NoError(tb, err)
	params.CircuitOrder = ModUpThenEncode
	params.ResidualParameters = residual
	require.Equal(tb, 32, params.EphemeralSecretWeight)
	require.Equal(tb, 13, params.BootstrappingParameters.LogN())
	require.Equal(tb, 8192, params.BootstrappingParameters.N())
	require.Equal(tb, 16, params.BootstrappingParameters.MaxLevel())
	require.Equal(tb, 45, residual.LogDefaultScale())
	require.Equal(tb, []uint64{
		36028797018652673, 549755731969, 549756026881, 549755486209, 549756174337,
		1152921504606830593, 1152921504606748673, 1152921504606994433,
		1152921504606683137, 1152921504606601217,
	}, params.BootstrappingParameters.Q()[:10], "Level 9 input and q9 divisor must match the retained helper fixture")
	qpJSON, err := json.Marshal(struct{ Q, P []uint64 }{params.BootstrappingParameters.Q(), params.BootstrappingParameters.P()})
	require.NoError(tb, err)
	qpSHA := sha256.Sum256(qpJSON)
	require.Equal(tb, batch027ExpectedQPSHA256, hex.EncodeToString(qpSHA[:]), "real E32 Q/P chain must match the pinned Primary profile")
	return params, residual
}

// The stream is sampled as three consecutive SplitMix64 words, most-significant
// word first, reduced modulo q0*q1*q2*q3, then centered into [-Q/2,Q/2].
// Generation order is component, coefficient, word; residues are emitted row 0..3.
func batch027SyntheticCiphertext(params Parameters, seed uint64) *rlwe.Ciphertext {
	const level, rows = batch027SourceLevel, batch027SourcePrefixRows
	ct := fastckks.NewCiphertext(params.BootstrappingParameters, 1, level)
	ct.Scale = rlwe.NewScale(uint64(1) << 45)
	ct.IsNTT, ct.IsMontgomery, ct.IsBatched = true, false, true
	ct.LogDimensions = ring.Dimensions{Cols: 12}
	q := params.BootstrappingParameters.Q()
	product := big.NewInt(1)
	for row := 0; row < rows; row++ {
		product.Mul(product, new(big.Int).SetUint64(q[row]))
	}
	half := new(big.Int).Rsh(new(big.Int).Set(product), 1)
	rng := batch027SplitMix64{state: seed}
	for component := range ct.Value {
		for coefficient := 0; coefficient < params.BootstrappingParameters.N(); coefficient++ {
			sample := new(big.Int)
			for word := 0; word < 3; word++ {
				sample.Lsh(sample, 64)
				sample.Or(sample, new(big.Int).SetUint64(rng.next()))
			}
			sample.Mod(sample, product)
			if sample.Cmp(half) > 0 {
				sample.Sub(sample, product)
			}
			for row := 0; row < rows; row++ {
				ct.Value[component].Coeffs[row][coefficient] = new(big.Int).Mod(new(big.Int).Set(sample), new(big.Int).SetUint64(q[row])).Uint64()
			}
		}
	}
	for component := range ct.Value {
		for row := 0; row < rows; row++ {
			params.BootstrappingParameters.RingQ().SubRings[row].NTT(ct.Value[component].Coeffs[row], ct.Value[component].Coeffs[row])
		}
	}
	return ct
}

func batch027QPHash(params Parameters) string {
	data, _ := json.Marshal(struct{ Q, P []uint64 }{params.BootstrappingParameters.Q(), params.BootstrappingParameters.P()})
	sum := sha256.Sum256(data)
	return hex.EncodeToString(sum[:])
}

// v1 hashes a domain/version separator, raw config SHA bytes, Q/P fingerprint
// plus Q/P arrays, LE metadata fields, canonical decimal scale, flags, and only
// authoritative c0/c1 q0..q3 NTT residues in component/row/coefficient order.
// Integers are uint64 LE; strings are UTF-8 with uint32-LE byte lengths.
func batch027FixtureDigest(ct *rlwe.Ciphertext, params Parameters) string {
	h := sha256.New()
	writeString := func(value string) {
		var size [4]byte
		binary.LittleEndian.PutUint32(size[:], uint32(len(value)))
		_, _ = h.Write(size[:])
		_, _ = h.Write([]byte(value))
	}
	writeU32 := func(value uint32) {
		var data [4]byte
		binary.LittleEndian.PutUint32(data[:], value)
		_, _ = h.Write(data[:])
	}
	_, _ = h.Write([]byte(batch027FixtureName + "\x00canonical-ntt-u64le-v1\x00"))
	configSHA, _ := hex.DecodeString(batch027ConfigSHA256)
	_, _ = h.Write(configSHA)
	writeString(batch027QPHash(params))
	q, p := params.BootstrappingParameters.Q(), params.BootstrappingParameters.P()
	writeU32(uint32(len(q)))
	for _, value := range q {
		var data [8]byte
		binary.LittleEndian.PutUint64(data[:], value)
		_, _ = h.Write(data[:])
	}
	writeU32(uint32(len(p)))
	for _, value := range p {
		var data [8]byte
		binary.LittleEndian.PutUint64(data[:], value)
		_, _ = h.Write(data[:])
	}
	for _, value := range []uint32{13, 8192, uint32(ct.Level()), uint32(ct.Degree()), 4, 8, 4, uint32(ct.LogDimensions.Rows), uint32(ct.LogDimensions.Cols)} {
		writeU32(value)
	}
	writeString(ct.Scale.BigInt().String())
	for _, flag := range []bool{ct.IsNTT, ct.IsMontgomery, ct.IsBatched} {
		if flag {
			_, _ = h.Write([]byte{1})
		} else {
			_, _ = h.Write([]byte{0})
		}
	}
	writeU32(uint32(len(ct.Value)))
	for component := range ct.Value {
		writeU32(4)
		for row := 0; row < batch027SourcePrefixRows; row++ {
			writeU32(uint32(len(ct.Value[component].Coeffs[row])))
			for _, coefficient := range ct.Value[component].Coeffs[row] {
				var data [8]byte
				binary.LittleEndian.PutUint64(data[:], coefficient)
				_, _ = h.Write(data[:])
			}
		}
	}
	return hex.EncodeToString(h.Sum(nil))
}

func batch027SameRows(a, b *rlwe.Ciphertext) bool {
	if len(a.Value) != len(b.Value) {
		return false
	}
	for component := range a.Value {
		for row := 0; row < batch027SourcePrefixRows; row++ {
			if len(a.Value[component].Coeffs[row]) != len(b.Value[component].Coeffs[row]) {
				return false
			}
			for coefficient := range a.Value[component].Coeffs[row] {
				if a.Value[component].Coeffs[row][coefficient] != b.Value[component].Coeffs[row][coefficient] {
					return false
				}
			}
		}
	}
	return true
}

func TestBatch027E32SyntheticFixtureReproducesAndMatchesGolden(t *testing.T) {
	params, _ := batch027E32Parameters(t)
	first := batch027SyntheticCiphertext(params, batch027FixtureSeed)
	second := batch027SyntheticCiphertext(params, batch027FixtureSeed)
	changedSeed := batch027SyntheticCiphertext(params, batch027FixtureSeed+1)
	require.True(t, batch027SameRows(first, second), "same-seed rebuilds must match every authoritative coefficient")
	require.Equal(t, first.Level(), second.Level())
	require.Equal(t, first.Degree(), second.Degree())
	require.Equal(t, first.Scale, second.Scale)
	require.Equal(t, first.IsNTT, second.IsNTT)
	require.Equal(t, first.IsMontgomery, second.IsMontgomery)
	require.Equal(t, first.IsBatched, second.IsBatched)
	require.Equal(t, first.LogDimensions, second.LogDimensions)
	firstSHA, secondSHA := batch027FixtureDigest(first, params), batch027FixtureDigest(second, params)
	require.Equal(t, firstSHA, secondSHA)
	t.Logf("fixture=%s seed=0x%016x qp_sha256=%s canonical_sha256=%s", batch027FixtureName, batch027FixtureSeed, batch027QPHash(params), firstSHA)
	require.Equal(t, batch027FixtureSHA256, firstSHA)
	require.NotEqual(t, firstSHA, batch027FixtureDigest(changedSeed, params), "different seed must change the fixture")
	require.Equal(t, batch027SourceLevel, first.Level())
	require.Equal(t, 1, first.Degree())
	require.True(t, first.Scale.Equal(rlwe.NewScale(uint64(1)<<45)))
	require.True(t, first.IsNTT)
	require.False(t, first.IsMontgomery)
	require.True(t, first.IsBatched)
	rows, err := fastckks.QPrefixWidth(first.Level())
	require.NoError(t, err)
	require.Equal(t, batch027SourcePrefixRows, rows)
	for component := range first.Value {
		nonzero := false
		for row := 0; row < batch027SourcePrefixRows; row++ {
			require.Len(t, first.Value[component].Coeffs[row], 8192)
			for _, coefficient := range first.Value[component].Coeffs[row] {
				nonzero = nonzero || coefficient != 0
			}
		}
		require.True(t, nonzero, "component %d must be nonzero", component)
		for row := batch027SourcePrefixRows; row < len(first.Value[component].Coeffs); row++ {
			require.Nil(t, first.Value[component].Coeffs[row], "dormant q%d must not have backing", row)
		}
	}
}

func batch027BigIntOracle(t testing.TB, params ckks.Parameters, input *rlwe.Ciphertext) [][][]uint64 {
	t.Helper()
	const sourceRows = batch027SourcePrefixRows
	q, ringQ := params.Q(), params.RingQ()
	qBig, inverses := make([]*big.Int, sourceRows), make([]*big.Int, sourceRows)
	product := big.NewInt(1)
	for row := 0; row < sourceRows; row++ {
		qBig[row] = new(big.Int).SetUint64(q[row])
		if row > 0 {
			productMod := new(big.Int).Mod(new(big.Int).Set(product), qBig[row])
			inverses[row] = new(big.Int).ModInverse(productMod, qBig[row])
			require.NotNil(t, inverses[row])
		}
		product.Mul(product, qBig[row])
	}
	divisor := new(big.Int).SetUint64(q[batch027SourceLevel])
	results := make([][][]uint64, len(input.Value))
	for component := range input.Value {
		coefficients := make([][]uint64, sourceRows)
		for row := 0; row < sourceRows; row++ {
			coefficients[row] = make([]uint64, params.N())
			ringQ.SubRings[row].INTT(input.Value[component].Coeffs[row], coefficients[row])
		}
		results[component] = make([][]uint64, batch027TargetRows)
		for row := 0; row < batch027TargetRows; row++ {
			results[component][row] = make([]uint64, params.N())
		}
		for coefficient := 0; coefficient < params.N(); coefficient++ {
			x, prefix := new(big.Int), big.NewInt(1)
			for row := 0; row < sourceRows; row++ {
				xMod := new(big.Int).Mod(new(big.Int).Set(x), qBig[row])
				delta := new(big.Int).Sub(new(big.Int).SetUint64(coefficients[row][coefficient]), xMod)
				delta.Mod(delta, qBig[row])
				digit := delta
				if row > 0 {
					digit.Mul(digit, inverses[row]).Mod(digit, qBig[row])
				}
				x.Add(x, new(big.Int).Mul(prefix, digit))
				prefix.Mul(prefix, qBig[row])
			}
			if x.Cmp(new(big.Int).Rsh(new(big.Int).Set(product), 1)) > 0 {
				x.Sub(x, product)
			}
			negative := x.Sign() < 0
			magnitude := new(big.Int).Abs(new(big.Int).Set(x))
			quotient, remainder := new(big.Int).QuoRem(magnitude, divisor, new(big.Int))
			if new(big.Int).Lsh(remainder, 1).Cmp(divisor) > 0 {
				quotient.Add(quotient, big.NewInt(1))
			}
			if negative {
				quotient.Neg(quotient)
			}
			twiceMagnitude := new(big.Int).Lsh(new(big.Int).Abs(new(big.Int).Set(quotient)), 1)
			require.Less(t, twiceMagnitude.Cmp(product), 0, "oracle output must fit Q0123 capacity")
			for row := 0; row < batch027TargetRows; row++ {
				results[component][row][coefficient] = new(big.Int).Mod(new(big.Int).Set(quotient), qBig[row]).Uint64()
			}
		}
		for row := 0; row < batch027TargetRows; row++ {
			ringQ.SubRings[row].NTT(results[component][row], results[component][row])
		}
	}
	return results
}

func TestBatch027E32PublicRescaleMatchesIndependentOracleInPlaceAndOutOfPlace(t *testing.T) {
	params, _ := batch027E32Parameters(t)
	eval, err := NewFastEvaluator(params)
	require.NoError(t, err)
	original := batch027SyntheticCiphertext(params, batch027FixtureSeed)
	want := batch027BigIntOracle(t, params.BootstrappingParameters, original)
	inputs := []*rlwe.Ciphertext{original.CopyNew(), original.CopyNew()}
	outputs := []*rlwe.Ciphertext{inputs[0], fastckks.NewCiphertext(params.BootstrappingParameters, 1, batch027TargetLevel)}
	outputs[1].Scale = rlwe.NewScale(1)
	for index := range inputs {
		require.NoError(t, eval.FastCKKS.Rescale(inputs[index], outputs[index]))
		require.Equal(t, batch027TargetLevel, outputs[index].Level())
		require.Equal(t, 1, outputs[index].Degree())
		wantScale := original.Scale.Div(rlwe.NewScale(params.BootstrappingParameters.Q()[batch027SourceLevel]))
		require.True(t, wantScale.Equal(outputs[index].Scale))
		require.True(t, outputs[index].IsNTT)
		require.False(t, outputs[index].IsMontgomery)
		for component := range outputs[index].Value {
			for row := 0; row < batch027TargetRows; row++ {
				require.Equal(t, want[component][row], outputs[index].Value[component].Coeffs[row], "in_place=%t component=%d q%d", index == 0, component, row)
			}
			for row := batch027TargetRows; row < len(outputs[index].Value[component].Coeffs); row++ {
				require.Nil(t, outputs[index].Value[component].Coeffs[row])
			}
		}
	}
}

func BenchmarkBatch027FastPublicRescaleE32SyntheticV1(b *testing.B) {
	params, _ := batch027E32Parameters(b)
	eval, err := NewFastEvaluator(params)
	require.NoError(b, err)
	original := batch027SyntheticCiphertext(params, batch027FixtureSeed)
	if got := batch027FixtureDigest(original, params); got != batch027FixtureSHA256 {
		b.Fatalf("fixture digest mismatch: got %s, want %s", got, batch027FixtureSHA256)
	}
	const batchSize = 64
	pool := make([]*rlwe.Ciphertext, batchSize)
	for i := range pool {
		pool[i] = original.CopyNew()
	}
	// Untimed evaluator-scratch warmup; each timed in-place call gets a fresh
	// copy of the immutable fixture, with resets excluded in 64-call batches.
	if err = eval.FastCKKS.Rescale(pool[0], pool[0]); err != nil {
		b.Fatal(err)
	}
	pool[0] = original.CopyNew()
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if i > 0 && i%batchSize == 0 {
			b.StopTimer()
			for j := range pool {
				pool[j] = original.CopyNew()
			}
			b.StartTimer()
		}
		ct := pool[i%batchSize]
		if err = eval.FastCKKS.Rescale(ct, ct); err != nil {
			b.Fatal(err)
		}
	}
}
