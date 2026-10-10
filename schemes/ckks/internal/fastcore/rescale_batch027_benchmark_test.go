package fastcore

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/ring"
)

const batch027CoreFixtureSeed = uint64(0x027e32c0ffee1234)

var (
	batch027CoreUint192Sink Uint192
	batch027CoreUint64Sink  uint64
)

type batch027CoreSplitMix64 struct{ state uint64 }

func (rng *batch027CoreSplitMix64) next() uint64 {
	// SplitMix64 v1; uint64 addition and multiplication wrap modulo 2^64.
	rng.state += 0x9e3779b97f4a7c15
	z := rng.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}

type batch027CoreSample struct {
	residues  [MaxQPrefixWidth]uint64
	magnitude Uint192
	negative  bool
}

func newBatch027CoreFixture(tb testing.TB) (*ring.Ring, *RescaleWorkspace, []batch027CoreSample) {
	tb.Helper()
	// Exact q0..q9 from the canonical Primary E32 profile. The companion
	// bootstrapping-package test asserts these values against its real Q/P SHA.
	q := []uint64{
		36028797018652673, 549755731969, 549756026881, 549755486209, 549756174337,
		1152921504606830593, 1152921504606748673, 1152921504606994433,
		1152921504606683137, 1152921504606601217,
	}
	ringQ, err := ring.NewRing(8192, q)
	require.NoError(tb, err)
	workspace := NewRescaleWorkspace(ringQ, 4)
	product := new(big.Int).SetUint64(workspace.modulus[3].Hi)
	product.Lsh(product, 64).Add(product, new(big.Int).SetUint64(workspace.modulus[3].Mid))
	product.Lsh(product, 64).Add(product, new(big.Int).SetUint64(workspace.modulus[3].Lo))
	half := new(big.Int).Rsh(new(big.Int).Set(product), 1)
	samples := make([]batch027CoreSample, 2*8192)
	rng := batch027CoreSplitMix64{state: batch027CoreFixtureSeed}
	// Flat sample order is component-major then coefficient; three generated
	// words concatenate most-significant first before modulo and centering.
	for i := range samples {
		value := new(big.Int)
		for word := 0; word < 3; word++ {
			value.Lsh(value, 64)
			value.Or(value, new(big.Int).SetUint64(rng.next()))
		}
		value.Mod(value, product)
		if value.Cmp(half) > 0 {
			value.Sub(value, product)
		}
		for row := 0; row < 4; row++ {
			samples[i].residues[row] = new(big.Int).Mod(new(big.Int).Set(value), new(big.Int).SetUint64(q[row])).Uint64()
		}
		reconstructed := reconstructQPrefix(4, samples[i].residues, workspace.q, workspace.inverse, workspace.modulus)
		samples[i].magnitude, samples[i].negative = centeredQPrefix(reconstructed, workspace.modulus[3], workspace.half[3])
	}
	return ringQ, workspace, samples
}

func TestBatch027CoreHelperFixtureUsesCanonicalQPrefix(t *testing.T) {
	ringQ, workspace, samples := newBatch027CoreFixture(t)
	require.Equal(t, 8192, ringQ.N())
	require.Equal(t, 9, ringQ.Level())
	require.Equal(t, 4, workspace.maxRows)
	require.Len(t, samples, 2*8192)
	for _, sample := range samples {
		rescaled := roundedMagnitude192(sample.magnitude, ringQ.SubRings[9].Modulus)
		if cmp192(rescaled, workspace.half[3]) >= 0 {
			t.Fatal("fixture rescale exceeded strict Q0123 capacity")
		}
	}
}

func BenchmarkBatch027FastCoreRescaleHelpers(b *testing.B) {
	ringQ, workspace, samples := newBatch027CoreFixture(b)
	q := []uint64{workspace.q[0], workspace.q[1], workspace.q[2], workspace.q[3]}
	divisor := ringQ.SubRings[9].Modulus
	b.Logf("fixture=batch027-e32-rescale-synthetic-v1 seed=0x%016x coefficient_samples=%d q9=%d; helper timings are non-additive and not full-Rescale percentages", batch027CoreFixtureSeed, len(samples), divisor)
	b.Run("reconstruct_q0123", func(b *testing.B) {
		var index int
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			batch027CoreUint192Sink = reconstructQPrefix(4, samples[index].residues, workspace.q, workspace.inverse, workspace.modulus)
			index++
			if index == len(samples) {
				index = 0
			}
		}
	})
	b.Run("rounded_magnitude_q9", func(b *testing.B) {
		var index int
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			batch027CoreUint192Sink = roundedMagnitude192(samples[index].magnitude, divisor)
			index++
			if index == len(samples) {
				index = 0
			}
		}
	})
	b.Run("signed_residue_q0123", func(b *testing.B) {
		var index int
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			sample := samples[index]
			batch027CoreUint64Sink = signedResidue192(sample.magnitude, sample.negative, q[index%len(q)])
			index++
			if index == len(samples) {
				index = 0
			}
		}
	})
}
