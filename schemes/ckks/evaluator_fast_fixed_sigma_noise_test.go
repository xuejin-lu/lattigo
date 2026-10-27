package ckks

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"math"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
)

func fixedSigmaTestRescaledOutput(t *testing.T) (Parameters, *Evaluator, *rlwe.Ciphertext) {
	return fixedSigmaTestRescaledOutputAtLevel(t, 4)
}

func fixedSigmaTestRescaledOutputAtLevel(t *testing.T, outputLevel int) (Parameters, *Evaluator, *rlwe.Ciphertext) {
	t.Helper()
	require.Contains(t, []int{3, 4}, outputLevel)
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39, 40, 39, 40, 60, 60}, LogP: []int{61},
		LogDefaultScale: 45, RingType: ring.Standard,
	})
	require.NoError(t, err)
	encoder := NewEncoder(params)
	sk := NewKeyGenerator(params).GenSecretKeyNew()
	encryptor := rlwe.NewEncryptor(params, sk)
	evaluator := NewEvaluator(params, nil)
	inputLevel := outputLevel + 1
	x := fixedSigmaTestEncrypt(t, params, encoder, encryptor, fastAddSubTestValues(2), rlwe.NewScale(math.Exp2(45)), inputLevel)
	y := fixedSigmaTestEncrypt(t, params, encoder, encryptor, fastAddSubTestValues(6), rlwe.NewScale(params.Q()[inputLevel]), inputLevel)
	product, err := evaluator.MulRelinNew(x, y)
	require.NoError(t, err)
	output := NewCiphertext(params, 1, outputLevel)
	require.NoError(t, evaluator.Rescale(product, output))
	require.Equal(t, outputLevel, output.Level())
	require.True(t, output.IsNTT)
	require.False(t, output.IsMontgomery)
	require.True(t, output.Scale.Equal(rlwe.NewScale(math.Exp2(45))))
	return params, evaluator, output
}

func fixedSigmaTestEncrypt(t *testing.T, params Parameters, encoder *Encoder, encryptor *rlwe.Encryptor, values []complex128, scale rlwe.Scale, level int) *rlwe.Ciphertext {
	t.Helper()
	plain := NewPlaintext(params, level)
	plain.Scale, plain.IsNTT = scale, true
	plain.LogDimensions = ring.Dimensions{Cols: 4}
	require.NoError(t, encoder.Encode(values, plain))
	ct, err := encryptor.EncryptNew(plain)
	require.NoError(t, err)
	return ct
}

func TestFixedSigmaCoefficientNoiseZeroIsStrictNoOp(t *testing.T) {
	_, evaluator, output := fixedSigmaTestRescaledOutput(t)
	before := snapshotCiphertextValue(output)
	evidence, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(output, 0, 17)
	require.NoError(t, err)
	require.True(t, evidence.NoOp)
	require.Empty(t, evidence.Phases)
	require.Zero(t, evidence.SampleCount)
	require.Empty(t, evidence.CoefficientSHA256)
	require.Equal(t, before, snapshotCiphertextValue(output))

	// The zero-sigma path deliberately does not inspect the receiver or target.
	evidence, err = (Evaluator{}).InjectFixedSigmaCoefficientNoiseDiagnostic(nil, 0, 99)
	require.NoError(t, err)
	require.True(t, evidence.NoOp)
}

func TestFixedSigmaCoefficientNoiseIsDeterministicAndMapsOnePolynomialToEveryRow(t *testing.T) {
	for _, level := range []int{4, 3} {
		t.Run(fmt.Sprintf("level_%d", level), func(t *testing.T) {
			params, evaluator, output := fixedSigmaTestRescaledOutputAtLevel(t, level)
			const sigma, seed = 3.75, int64(20261022)
			first, second := output.CopyNew(), output.CopyNew()
			before := first.CopyNew()
			c1Before := snapshotCiphertextValueComponent(first, 1)
			dormantBefore := dormantQ4Backing(first)
			evidenceFirst, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(first, sigma, seed)
			require.NoError(t, err)
			evidenceSecond, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(second, sigma, seed)
			require.NoError(t, err)
			firstEvidenceWithoutTiming, secondEvidenceWithoutTiming := evidenceFirst, evidenceSecond
			firstEvidenceWithoutTiming.Phases, secondEvidenceWithoutTiming.Phases = nil, nil
			require.Equal(t, firstEvidenceWithoutTiming, secondEvidenceWithoutTiming)
			require.Equal(t, []string{
				"sampling_hash_moments", "coefficient_to_qprefix_ntt", "candidate_c0_add",
				"qprefix_centered_capacity_check", "commit_c0",
			}, fixedSigmaCoefficientNoisePhaseNames(evidenceFirst.Phases))
			for _, phase := range evidenceFirst.Phases {
				require.GreaterOrEqual(t, phase.ElapsedNS, int64(0))
				require.True(t, phase.MemoryStatsCollected, "historical helper keeps its existing allocation diagnostics")
			}
			require.Equal(t, first.Value[0], second.Value[0])
			require.Equal(t, evidenceFirst.SampleCount, params.N())
			require.Equal(t, 4, evidenceFirst.Capacity.Rows)
			require.Equal(t, params.N(), evidenceFirst.Capacity.CheckedCoefficients)
			require.True(t, evidenceFirst.Capacity.StrictFit)
			require.NotEmpty(t, evidenceFirst.Capacity.PrefixProduct)
			require.NotEmpty(t, evidenceFirst.CoefficientSHA256)
			decodedHash, err := hex.DecodeString(evidenceFirst.CoefficientSHA256)
			require.NoError(t, err)
			require.Len(t, decodedHash, 32)

			// Only c0 changes. Ciphertext metadata, c1, and the dormant q4 backing stay
			// identical; the test never inspects dormant row contents.
			require.Equal(t, before.Degree(), first.Degree())
			require.Equal(t, before.Level(), first.Level())
			require.Equal(t, *before.MetaData, *first.MetaData)
			require.True(t, first.Scale.Equal(rlwe.NewScale(math.Exp2(45))))
			require.Equal(t, before.IsNTT, first.IsNTT)
			require.Equal(t, before.IsMontgomery, first.IsMontgomery)
			require.Equal(t, c1Before, snapshotCiphertextValueComponent(first, 1))
			require.Equal(t, dormantBefore, dormantQ4Backing(first))
			require.Empty(t, dormantBefore)

			// Recreate the specified rounded Gaussian sequence, then verify the
			// inverse-transformed delta in every maintained residue row equals that
			// same signed integer polynomial modulo the row's q_i.
			rng := rand.New(rand.NewSource(seed))
			expected := make([]int64, params.N())
			for coefficient := range expected {
				expected[coefficient] = int64(math.Round(rng.NormFloat64() * sigma))
			}
			for row := 0; row < 4; row++ {
				subring := params.RingQ().SubRings[row]
				q := subring.Modulus
				deltaNTT := make([]uint64, params.N())
				for coefficient := range deltaNTT {
					after, prior := first.Value[0].Coeffs[row][coefficient], before.Value[0].Coeffs[row][coefficient]
					if after >= prior {
						deltaNTT[coefficient] = after - prior
					} else {
						deltaNTT[coefficient] = q - (prior - after)
					}
				}
				coefficientDelta := make([]uint64, params.N())
				subring.INTT(deltaNTT, coefficientDelta)
				for coefficient, got := range coefficientDelta {
					require.Equal(t, signedCoefficientResidue(expected[coefficient], q), got,
						"q%d coefficient %d must encode the same sampled signed integer", row, coefficient)
				}
			}
		})
	}
}

func fixedSigmaCoefficientNoisePhaseNames(phases []FixedSigmaCoefficientNoisePhase) []string {
	names := make([]string, len(phases))
	for index, phase := range phases {
		names[index] = phase.Name
	}
	return names
}

func TestFixedSigmaCoefficientNoiseRejectsNonzeroC1WithoutMutation(t *testing.T) {
	_, evaluator, output := fixedSigmaTestRescaledOutput(t)
	output.Value[1].Coeffs[1][0] = 1
	before := snapshotCiphertextValue(output)
	metadata := *output.MetaData
	_, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(output, 3.75, 9)
	require.ErrorContains(t, err, "c1=0")
	require.Equal(t, before, snapshotCiphertextValue(output))
	require.Equal(t, metadata, *output.MetaData)
}

func TestFixedSigmaCoefficientNoiseRejectsLevelsOutsideBoundedProfile(t *testing.T) {
	params, evaluator, _ := fixedSigmaTestRescaledOutput(t)
	for _, level := range []int{2, 5} {
		t.Run(fmt.Sprintf("level_%d", level), func(t *testing.T) {
			unsupported := NewCiphertext(params, 1, level)
			unsupported.Scale = rlwe.NewScale(math.Exp2(45))
			unsupported.IsNTT = true
			_, err := evaluator.InjectFixedSigmaCoefficientNoiseDiagnostic(unsupported, 3.75, 9)
			require.ErrorContains(t, err, "Level4 or Level3")
		})
	}
}

func TestFixedSigmaBootstrapMulNoiseSupportsActiveLevelsAndScales(t *testing.T) {
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39, 40, 39, 40, 60, 60}, LogP: []int{61},
		LogDefaultScale: 45, RingType: ring.Standard,
	})
	require.NoError(t, err)
	evaluator := NewEvaluator(params, nil)
	for _, level := range []int{0, 1, 2, 3, 4} {
		for _, scale := range []float64{math.Exp2(30), math.Exp2(45), math.Exp2(60)} {
			t.Run(fmt.Sprintf("level_%d_scale_%.0f", level, math.Log2(scale)), func(t *testing.T) {
				ct := NewCiphertext(params, 1, level)
				ct.Scale, ct.IsNTT = rlwe.NewScale(scale), true
				before := ct.CopyNew()
				evidence, err := evaluator.InjectFixedSigmaBootstrapMulNoiseDiagnostic(ct, 4.004280012097432, int64(level)+17)
				require.NoError(t, err)
				require.Equal(t, level, evidence.Level)
				require.Equal(t, ct.Scale.Value.Text('g', -1), evidence.Scale)
				require.Equal(t, min(level+1, 4), evidence.Capacity.Rows)
				require.True(t, evidence.Capacity.StrictFit)
				for _, phase := range evidence.Phases {
					require.False(t, phase.MemoryStatsCollected, "Bootstrap hot-path helper must not call runtime.ReadMemStats per phase")
					require.Zero(t, phase.AllocBytes)
					require.Zero(t, phase.Allocs)
				}
				require.Equal(t, before.Level(), ct.Level())
				require.Equal(t, *before.MetaData, *ct.MetaData)
				require.True(t, before.Scale.Equal(ct.Scale))
				require.Equal(t, before.Value[1], ct.Value[1])
				require.NotEqual(t, before.Value[0], ct.Value[0])
			})
		}
	}
}

func TestFixedSigmaBootstrapMulNoiseMontgomeryLevel12Scale60PreservesSignedCoefficients(t *testing.T) {
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: 13,
		LogQ: []int{55, 39, 40, 39, 40, 60, 60, 60, 60, 60, 60, 60, 60, 56, 57, 56, 56},
		LogP: []int{61, 61, 61, 61, 61}, LogDefaultScale: 45, RingType: ring.Standard,
	})
	require.NoError(t, err)
	evaluator := NewEvaluator(params, nil)
	const (
		level = 12
		sigma = 4.004280012097432
		seed  = int64(20261027)
	)
	scale := rlwe.NewScale(math.Exp2(60))
	ct := NewCiphertext(params, 1, level)
	ct.Scale, ct.IsNTT, ct.IsMontgomery = scale, true, true
	ringQ := params.RingQ()
	rows := 4

	// Encode the same small signed coefficient polynomial into each active Q
	// row, then put it in the genuine NTT+Montgomery representation used by the
	// Fast Chebyshev power path.
	for row := 0; row < rows; row++ {
		for i := range ct.Value[0].Coeffs[row] {
			coefficient := int64(i%7 - 3)
			ct.Value[0].Coeffs[row][i] = signedCoefficientResidue(coefficient, ringQ.SubRings[row].Modulus)
		}
		ringQ.SubRings[row].NTT(ct.Value[0].Coeffs[row], ct.Value[0].Coeffs[row])
		ringQ.SubRings[row].MForm(ct.Value[0].Coeffs[row], ct.Value[0].Coeffs[row])
	}
	before := ct.CopyNew()
	beforeMeta := *ct.MetaData

	evidence, err := evaluator.InjectFixedSigmaBootstrapMulNoiseDiagnostic(ct, sigma, seed)
	require.NoError(t, err)
	require.Equal(t, level, evidence.Level)
	require.True(t, evidence.Scale == ct.Scale.Value.Text('g', -1))
	require.True(t, evidence.Capacity.StrictFit)
	require.Equal(t, rows, evidence.Capacity.Rows)
	require.Equal(t, level, ct.Level())
	require.Equal(t, beforeMeta, *ct.MetaData)
	require.True(t, ct.Scale.Equal(scale))
	require.True(t, ct.IsNTT)
	require.True(t, ct.IsMontgomery)
	require.Equal(t, before.Value[1], ct.Value[1], "c1 must be preserved")
	for row := rows; row <= level; row++ {
		require.Equal(t, before.Value[0].Coeffs[row], ct.Value[0].Coeffs[row], "dormant q%d must be preserved", row)
	}

	// Convert copies back to signed-coefficient residues and verify every active
	// row contains exactly the original coefficient plus the sampled Gaussian.
	gotCoefficients := ct.Value[0].CopyNew()
	beforeCoefficients := before.Value[0].CopyNew()
	for row := 0; row < rows; row++ {
		subring := ringQ.SubRings[row]
		subring.INTT(gotCoefficients.Coeffs[row], gotCoefficients.Coeffs[row])
		subring.IMForm(gotCoefficients.Coeffs[row], gotCoefficients.Coeffs[row])
		subring.INTT(beforeCoefficients.Coeffs[row], beforeCoefficients.Coeffs[row])
		subring.IMForm(beforeCoefficients.Coeffs[row], beforeCoefficients.Coeffs[row])
	}
	rng := rand.New(rand.NewSource(seed))
	for coefficient := 0; coefficient < params.N(); coefficient++ {
		noise := int64(math.Round(rng.NormFloat64() * sigma))
		for row := 0; row < rows; row++ {
			modulus := ringQ.SubRings[row].Modulus
			want := (beforeCoefficients.Coeffs[row][coefficient] + signedCoefficientResidue(noise, modulus)) % modulus
			require.Equalf(t, want, gotCoefficients.Coeffs[row][coefficient], "q%d coefficient %d", row, coefficient)
		}
	}
	require.Contains(t, evidence.DispatchEvidence, "MForm-mapped")
}

func fixedSigmaBootstrapMulNoiseLogN16Parameters(t *testing.T) Parameters {
	return fixedSigmaBootstrapMulNoiseParametersForLogN(t, 16)
}

func fixedSigmaBootstrapMulNoiseParametersForLogN(t *testing.T, logN int) Parameters {
	t.Helper()
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: logN,
		Q: []uint64{
			36028797019488257, 549754109953, 549753978881, 549753716737, 549753192449,
			1152921504606584833, 1152921504614055937, 1152921504598720513, 1152921504615628801,
			1152921504616808449, 1152921504597016577, 1152921504595968001, 1152921504618381313,
			72057594038321153, 72057594036879361, 72057594035306497, 72057594040680449,
		},
		P:               []uint64{2305843009211596801, 2305843009210023937, 2305843009218281473, 2305843009208713217, 2305843009218936833},
		LogDefaultScale: 45, RingType: ring.Standard,
	})
	require.NoError(t, err)
	require.Equal(t, 1<<logN, params.N())
	require.Len(t, params.Q(), 17)
	require.Len(t, params.P(), 5)
	return params
}

func TestFixedSigmaBootstrapMulNoiseLogN16MontgomeryLevel12Scale60PreservesSignedCoefficients(t *testing.T) {
	params := fixedSigmaBootstrapMulNoiseLogN16Parameters(t)
	evaluator := NewEvaluator(params, nil)
	const (
		level = 12
		sigma = 4.004280012097432
		seed  = int64(20261027)
	)
	scale := rlwe.NewScale(math.Exp2(60))
	ct := NewCiphertext(params, 1, level)
	ct.Scale, ct.IsNTT, ct.IsMontgomery = scale, true, true
	ringQ := params.RingQ()
	const rows = 4
	for row := 0; row < rows; row++ {
		for i := range ct.Value[0].Coeffs[row] {
			ct.Value[0].Coeffs[row][i] = signedCoefficientResidue(int64(i%7-3), ringQ.SubRings[row].Modulus)
		}
		ringQ.SubRings[row].NTT(ct.Value[0].Coeffs[row], ct.Value[0].Coeffs[row])
		ringQ.SubRings[row].MForm(ct.Value[0].Coeffs[row], ct.Value[0].Coeffs[row])
	}
	before := ct.CopyNew()
	beforeMeta := *ct.MetaData

	evidence, err := evaluator.InjectFixedSigmaBootstrapMulNoiseDiagnostic(ct, sigma, seed)
	require.NoError(t, err)
	require.Equal(t, level, evidence.Level)
	require.Equal(t, scale.Value.Text('g', -1), evidence.Scale)
	require.True(t, evidence.Capacity.StrictFit)
	require.Equal(t, rows, evidence.Capacity.Rows)
	require.Equal(t, params.N(), evidence.Capacity.CheckedCoefficients)
	require.Equal(t, beforeMeta, *ct.MetaData)
	require.True(t, ct.Scale.Equal(scale))
	require.True(t, ct.IsNTT)
	require.True(t, ct.IsMontgomery)
	require.Equal(t, before.Value[1], ct.Value[1])
	for row := rows; row <= level; row++ {
		require.Equal(t, before.Value[0].Coeffs[row], ct.Value[0].Coeffs[row], "dormant q%d must be unchanged", row)
	}

	got := ct.Value[0].CopyNew()
	prior := before.Value[0].CopyNew()
	for row := 0; row < rows; row++ {
		ringQ.SubRings[row].INTT(got.Coeffs[row], got.Coeffs[row])
		ringQ.SubRings[row].IMForm(got.Coeffs[row], got.Coeffs[row])
		ringQ.SubRings[row].INTT(prior.Coeffs[row], prior.Coeffs[row])
		ringQ.SubRings[row].IMForm(prior.Coeffs[row], prior.Coeffs[row])
	}
	rng := rand.New(rand.NewSource(seed))
	for coefficient := 0; coefficient < params.N(); coefficient++ {
		noise := int64(math.Round(rng.NormFloat64() * sigma))
		for row := 0; row < rows; row++ {
			modulus := ringQ.SubRings[row].Modulus
			want := (prior.Coeffs[row][coefficient] + signedCoefficientResidue(noise, modulus)) % modulus
			require.Equalf(t, want, got.Coeffs[row][coefficient], "q%d coefficient %d", row, coefficient)
		}
	}
}

func TestFixedSigmaBootstrapMulNoiseLogN16CapacityFailureDoesNotMutateTarget(t *testing.T) {
	params := fixedSigmaBootstrapMulNoiseLogN16Parameters(t)
	evaluator := NewEvaluator(params, nil)
	ct := NewCiphertext(params, 1, 0)
	ct.Scale, ct.IsNTT = rlwe.NewScale(math.Exp2(45)), true
	q := params.RingQ().SubRings[0].Modulus
	seed := int64(-1)
	for candidateSeed := int64(0); candidateSeed < 1000; candidateSeed++ {
		firstNoise := int64(math.Round(rand.New(rand.NewSource(candidateSeed)).NormFloat64() * 4.004280012097432))
		if firstNoise > 1 {
			seed = candidateSeed
			break
		}
	}
	require.GreaterOrEqual(t, seed, int64(0))
	coefficient := make([]uint64, params.N())
	coefficient[0] = q/2 - 1
	params.RingQ().SubRings[0].NTT(coefficient, ct.Value[0].Coeffs[0])
	before := snapshotCiphertextValue(ct)
	metadata := *ct.MetaData
	_, err := evaluator.InjectFixedSigmaBootstrapMulNoiseDiagnostic(ct, 4.004280012097432, seed)
	require.ErrorContains(t, err, "exceeds or cannot prove Q-prefix capacity")
	require.Equal(t, before, snapshotCiphertextValue(ct))
	require.Equal(t, metadata, *ct.MetaData)
}

func TestFixedSigmaBootstrapMulNoiseSupportsLogN13Through16MontgomeryEncodedInputs(t *testing.T) {
	const (
		level = 12
		sigma = 4.004280012097432
	)
	for _, logN := range []int{13, 14, 15, 16} {
		t.Run(fmt.Sprintf("logn%d", logN), func(t *testing.T) {
			params := fixedSigmaBootstrapMulNoiseParametersForLogN(t, logN)
			encoder := NewEncoder(params)
			plain := NewPlaintext(params, level)
			plain.Scale = rlwe.NewScale(math.Exp2(60))
			plain.IsNTT = true
			values := make([]complex128, params.MaxDimensions().Cols)
			values[0], values[1] = complex(0.125, -0.25), complex(-0.0625, 0.03125)
			require.NoError(t, encoder.Encode(values, plain))

			ct := NewCiphertext(params, 1, level)
			ct.Value[0] = *plain.Value.CopyNew()
			ct.Scale = plain.Scale
			ct.IsNTT, ct.IsBatched = plain.IsNTT, plain.IsBatched
			ct.LogDimensions, ct.IsBitReversed = plain.LogDimensions, plain.IsBitReversed
			ct.IsMontgomery = true
			for row := 0; row < 4; row++ {
				params.RingQ().SubRings[row].MForm(ct.Value[0].Coeffs[row], ct.Value[0].Coeffs[row])
			}
			before := ct.CopyNew()
			beforeMeta := *ct.MetaData
			seed := int64(20261027 + logN)
			evaluator := NewEvaluator(params, nil)
			evidence, err := evaluator.InjectFixedSigmaBootstrapMulNoiseDiagnostic(ct, sigma, seed)
			require.NoError(t, err)
			require.Equal(t, params.N(), evidence.SampleCount)
			require.Equal(t, seed, evidence.Seed)
			require.Equal(t, 4, evidence.Capacity.Rows)
			require.Equal(t, params.N(), evidence.Capacity.CheckedCoefficients)
			require.True(t, evidence.Capacity.StrictFit)
			require.Len(t, mustDecodeSHA256(t, evidence.CoefficientSHA256), sha256.Size)
			require.Equal(t, beforeMeta, *ct.MetaData)
			require.Equal(t, before.Value[1], ct.Value[1], "c1 must be unchanged")
			require.Equal(t, level, ct.Level())
			require.True(t, ct.Scale.Equal(rlwe.NewScale(math.Exp2(60))))
			require.True(t, ct.IsNTT && ct.IsMontgomery)
			for row := 4; row <= level; row++ {
				require.Equal(t, before.Value[0].Coeffs[row], ct.Value[0].Coeffs[row], "dormant q%d must be unchanged", row)
			}

			// The evidence digest must be exactly the deterministic signed Gaussian
			// sequence of length params.N(), and each active Montgomery/NTT row must
			// decode to the original encoded coefficient plus that same sequence.
			rng := rand.New(rand.NewSource(seed))
			expected := make([]int64, params.N())
			hasher := sha256.New()
			var encoded [8]byte
			for index := range expected {
				expected[index] = int64(math.Round(rng.NormFloat64() * sigma))
				binary.LittleEndian.PutUint64(encoded[:], uint64(expected[index]))
				_, _ = hasher.Write(encoded[:])
			}
			require.Equal(t, hex.EncodeToString(hasher.Sum(nil)), evidence.CoefficientSHA256)
			got, prior := ct.Value[0].CopyNew(), before.Value[0].CopyNew()
			for row := 0; row < evidence.Capacity.Rows; row++ {
				subring := params.RingQ().SubRings[row]
				subring.INTT(got.Coeffs[row], got.Coeffs[row])
				subring.IMForm(got.Coeffs[row], got.Coeffs[row])
				subring.INTT(prior.Coeffs[row], prior.Coeffs[row])
				subring.IMForm(prior.Coeffs[row], prior.Coeffs[row])
				modulus := subring.Modulus
				for coefficient, noise := range expected {
					want := (prior.Coeffs[row][coefficient] + signedCoefficientResidue(noise, modulus)) % modulus
					require.Equalf(t, want, got.Coeffs[row][coefficient], "q%d coefficient %d", row, coefficient)
				}
			}

			// A near-boundary Level-0 coefficient must still fail the strict
			// capacity proof transactionally at every ring dimension.
			failure := NewCiphertext(params, 1, 0)
			failure.Scale = rlwe.NewScale(math.Exp2(60))
			failure.IsNTT, failure.IsMontgomery = true, true
			q0 := params.RingQ().SubRings[0]
			failure.Value[0].Coeffs[0][0] = q0.Modulus/2 - 1
			q0.NTT(failure.Value[0].Coeffs[0], failure.Value[0].Coeffs[0])
			q0.MForm(failure.Value[0].Coeffs[0], failure.Value[0].Coeffs[0])
			failureSeed := fixedSigmaCapacityFailureSeed(sigma)
			require.GreaterOrEqual(t, failureSeed, int64(0))
			failureBefore := snapshotCiphertextValue(failure)
			_, err = evaluator.InjectFixedSigmaBootstrapMulNoiseDiagnostic(failure, sigma, failureSeed)
			require.ErrorContains(t, err, "exceeds or cannot prove Q-prefix capacity")
			require.Equal(t, failureBefore, snapshotCiphertextValue(failure))
		})
	}
}

func mustDecodeSHA256(t *testing.T, value string) []byte {
	t.Helper()
	digest, err := hex.DecodeString(value)
	require.NoError(t, err)
	return digest
}

func fixedSigmaCapacityFailureSeed(sigma float64) int64 {
	for seed := int64(0); seed < 1000; seed++ {
		if int64(math.Round(rand.New(rand.NewSource(seed)).NormFloat64()*sigma)) > 1 {
			return seed
		}
	}
	return -1
}

func TestFixedSigmaBootstrapMulNoiseCapacityFailureDoesNotMutateTarget(t *testing.T) {
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39, 40, 39, 40, 60, 60}, LogP: []int{61},
		LogDefaultScale: 45, RingType: ring.Standard,
	})
	require.NoError(t, err)
	evaluator := NewEvaluator(params, nil)
	ct := NewCiphertext(params, 1, 0)
	ct.Scale, ct.IsNTT = rlwe.NewScale(math.Exp2(37)), true
	q := params.RingQ().SubRings[0].Modulus
	var firstNoise int64
	var seed int64
	for candidateSeed := int64(0); candidateSeed < 1000; candidateSeed++ {
		firstNoise = int64(math.Round(rand.New(rand.NewSource(candidateSeed)).NormFloat64() * 4.004280012097432))
		if firstNoise > 1 {
			seed = candidateSeed
			break
		}
	}
	require.Greater(t, firstNoise, int64(1))
	coefficient := make([]uint64, params.N())
	coefficient[0] = q/2 - 1
	params.RingQ().SubRings[0].NTT(coefficient, ct.Value[0].Coeffs[0])
	before := snapshotCiphertextValue(ct)
	meta := *ct.MetaData
	_, err = evaluator.InjectFixedSigmaBootstrapMulNoiseDiagnostic(ct, 4.004280012097432, seed)
	require.ErrorContains(t, err, "exceeds or cannot prove Q-prefix capacity")
	require.Equal(t, before, snapshotCiphertextValue(ct))
	require.Equal(t, meta, *ct.MetaData)
}

func TestFixedSigmaBootstrapMulNoiseMontgomeryCapacityFailureDoesNotMutateTarget(t *testing.T) {
	params, err := NewParametersFromLiteral(ParametersLiteral{
		LogN: 13, LogQ: []int{55, 39, 40, 39, 40, 60, 60}, LogP: []int{61},
		LogDefaultScale: 45, RingType: ring.Standard,
	})
	require.NoError(t, err)
	evaluator := NewEvaluator(params, nil)
	ct := NewCiphertext(params, 1, 0)
	ct.Scale, ct.IsNTT, ct.IsMontgomery = rlwe.NewScale(math.Exp2(60)), true, true
	q := params.RingQ().SubRings[0].Modulus
	var firstNoise int64
	var seed int64
	for candidateSeed := int64(0); candidateSeed < 1000; candidateSeed++ {
		firstNoise = int64(math.Round(rand.New(rand.NewSource(candidateSeed)).NormFloat64() * 4.004280012097432))
		if firstNoise > 1 {
			seed = candidateSeed
			break
		}
	}
	require.Greater(t, firstNoise, int64(1))
	coefficient := make([]uint64, params.N())
	coefficient[0] = q/2 - 1
	subring := params.RingQ().SubRings[0]
	subring.NTT(coefficient, ct.Value[0].Coeffs[0])
	subring.MForm(ct.Value[0].Coeffs[0], ct.Value[0].Coeffs[0])
	before := snapshotCiphertextValue(ct)
	meta := *ct.MetaData
	_, err = evaluator.InjectFixedSigmaBootstrapMulNoiseDiagnostic(ct, 4.004280012097432, seed)
	require.ErrorContains(t, err, "exceeds or cannot prove Q-prefix capacity")
	require.Equal(t, before, snapshotCiphertextValue(ct))
	require.Equal(t, meta, *ct.MetaData)
}

func snapshotCiphertextValue(ct *rlwe.Ciphertext) [][][]uint64 {
	values := make([][][]uint64, len(ct.Value))
	for component := range ct.Value {
		values[component] = snapshotCiphertextValueComponent(ct, component)
	}
	return values
}

func snapshotCiphertextValueComponent(ct *rlwe.Ciphertext, component int) [][]uint64 {
	poly := ct.Value[component]
	if poly.Coeffs == nil {
		return nil
	}
	rows := make([][]uint64, len(poly.Coeffs))
	for row, values := range poly.Coeffs {
		if values != nil {
			rows[row] = append([]uint64{}, values...)
		}
	}
	return rows
}

func rowBackingLengths(values []uint64) []int {
	if len(values) == 0 {
		return nil
	}
	return []int{len(values)}
}

func dormantQ4Backing(ct *rlwe.Ciphertext) []int {
	if len(ct.Value[0].Coeffs) <= 4 {
		return nil
	}
	return rowBackingLengths(ct.Value[0].Coeffs[4])
}
