package ckks

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"math"
	"math/rand"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

// FixedSigmaCoefficientNoiseEvidence describes one explicit diagnostic-only
// fixed-sigma injection. Coefficients are rounded signed Gaussian integers;
// the same signed polynomial is reduced into every maintained Q row.
type FixedSigmaCoefficientNoiseEvidence struct {
	Sigma               float64                          `json:"sigma_coefficient"`
	Seed                int64                            `json:"sampling_seed"`
	NoOp                bool                             `json:"no_op"`
	SampleCount         int                              `json:"coefficient_sample_count"`
	CoefficientMean     float64                          `json:"coefficient_mean"`
	CoefficientVariance float64                          `json:"coefficient_variance"`
	CoefficientRMS      float64                          `json:"coefficient_rms"`
	CoefficientMaxAbs   uint64                           `json:"coefficient_max_abs"`
	CoefficientSHA256   string                           `json:"coefficient_sha256"`
	SamplingAlgorithm   string                           `json:"sampling_algorithm"`
	Capacity            fastcore.QPrefixCapacityEvidence `json:"qprefix_capacity"`
	DispatchEvidence    string                           `json:"dispatch_evidence"`
}

// InjectFixedSigmaCoefficientNoiseDiagnostic explicitly adds one sampled
// coefficient-domain polynomial to a Fast zero-secret ciphertext's c0. It is
// diagnostic-only, supports the isolated LogN13 task boundary (Level 4,
// four-row Q-prefix), and is never called by ordinary CKKS operations.
func (eval Evaluator) InjectFixedSigmaCoefficientNoiseDiagnostic(ct *rlwe.Ciphertext, sigma float64, seed int64) (FixedSigmaCoefficientNoiseEvidence, error) {
	evidence := FixedSigmaCoefficientNoiseEvidence{Sigma: sigma, Seed: seed}
	if math.IsNaN(sigma) || math.IsInf(sigma, 0) || sigma < 0 {
		return evidence, fmt.Errorf("fixed-sigma diagnostic requires finite sigma >= 0")
	}
	// A zero sigma is a strict no-op: it does not inspect or mutate ct and does
	// not instantiate or advance a random source.
	if sigma == 0 {
		evidence.NoOp = true
		evidence.SamplingAlgorithm = "NOT_RUN_ZERO_SIGMA"
		evidence.DispatchEvidence = "sigma=0; ciphertext and sampler untouched"
		return evidence, nil
	}
	if eval.Encoder == nil || eval.Evaluator == nil {
		return evidence, fmt.Errorf("fixed-sigma diagnostic requires an initialized CKKS evaluator")
	}
	if !eval.usesFastCKKSZeroSecretSimulation() {
		return evidence, fmt.Errorf("fixed-sigma diagnostic is available only for the Fast zero-secret CKKS capability")
	}
	params := eval.GetParameters()
	if _, err := validateFastCKKSZeroSecretMulShape("fixed-sigma target", ct, params); err != nil {
		return evidence, fmt.Errorf("fixed-sigma diagnostic target: %w", err)
	}
	if params.LogN() != 13 || params.N() != 1<<13 || ct.Level() != 4 || !ct.Scale.Equal(rlwe.NewScale(math.Exp2(45))) ||
		!ct.IsNTT || ct.IsMontgomery || !eval.GetRLWEParameters().NTTFlag() {
		return evidence, fmt.Errorf("fixed-sigma diagnostic requires the LogN13, Level4, Scale2^45 NTT non-Montgomery output profile")
	}
	rows, err := fastcore.QPrefixWidth(ct.Level())
	if err != nil || rows != 4 {
		return evidence, fmt.Errorf("fixed-sigma diagnostic requires exactly four maintained Q-prefix rows: rows=%d err=%v", rows, err)
	}
	ringQ := params.RingQ()
	if err := fastcore.ValidatePrefixRows(ringQ, ct.Level(), rows, ct.Value...); err != nil {
		return evidence, fmt.Errorf("fixed-sigma diagnostic Q-prefix backing: %w", err)
	}
	for row := 0; row < rows; row++ {
		modulus := ringQ.SubRings[row].Modulus
		for coefficient, value := range ct.Value[0].Coeffs[row] {
			if value >= modulus {
				return evidence, fmt.Errorf("fixed-sigma diagnostic c0 q%d coefficient %d is not a canonical residue", row, coefficient)
			}
		}
		for coefficient, value := range ct.Value[1].Coeffs[row] {
			if value != 0 {
				return evidence, fmt.Errorf("fixed-sigma diagnostic requires c1=0; q%d coefficient %d is nonzero", row, coefficient)
			}
		}
	}

	n := params.N()
	evidence.SamplingAlgorithm = "Go math/rand.NewSource(seed).NormFloat64(), rounded to nearest integer with math.Round (ties away from zero)"
	coefficients := make([]int64, n)
	rng := rand.New(rand.NewSource(seed))
	hasher := sha256.New()
	var mean, m2, sumSquares float64
	var maxAbs uint64
	var encoded [8]byte
	for index := range coefficients {
		sampled := math.Round(rng.NormFloat64() * sigma)
		if math.IsNaN(sampled) || math.IsInf(sampled, 0) || sampled >= math.Exp2(63) || sampled < -math.Exp2(63) {
			return evidence, fmt.Errorf("fixed-sigma diagnostic sample %d cannot be represented as int64", index)
		}
		value := int64(sampled)
		coefficients[index] = value
		binary.LittleEndian.PutUint64(encoded[:], uint64(value))
		_, _ = hasher.Write(encoded[:])
		abs := uint64(value)
		if value < 0 {
			abs = uint64(-(value + 1)) + 1
		}
		if abs > maxAbs {
			maxAbs = abs
		}
		observation := float64(value)
		delta := observation - mean
		mean += delta / float64(index+1)
		m2 += delta * (observation - mean)
		sumSquares += observation * observation
	}
	variance := m2 / float64(n)
	evidence.SampleCount = n
	evidence.CoefficientMean = mean
	evidence.CoefficientVariance = variance
	evidence.CoefficientRMS = math.Sqrt(sumSquares / float64(n))
	evidence.CoefficientMaxAbs = maxAbs
	evidence.CoefficientSHA256 = hex.EncodeToString(hasher.Sum(nil))

	noiseNTT := ring.NewPoly(n, rows-1)
	for row := 0; row < rows; row++ {
		modulus := ringQ.SubRings[row].Modulus
		for coefficient, value := range coefficients {
			noiseNTT.Coeffs[row][coefficient] = signedCoefficientResidue(value, modulus)
		}
		ringQ.SubRings[row].NTT(noiseNTT.Coeffs[row], noiseNTT.Coeffs[row])
	}

	candidate := ct.Value[0].CopyNew()
	for row := 0; row < rows; row++ {
		ringQ.SubRings[row].Add(candidate.Coeffs[row], noiseNTT.Coeffs[row], candidate.Coeffs[row])
	}
	candidateCT := *ct
	candidateCT.Value = append([]ring.Poly(nil), ct.Value...)
	candidateCT.Value[0] = *candidate
	capacityChecker, ok := eval.rescaleCore.(interface {
		CheckQPrefixAdditionCapacity(*ring.Ring, *rlwe.Ciphertext, *rlwe.Ciphertext, []int64, int) (fastcore.QPrefixCapacityEvidence, error)
	})
	if !ok {
		return evidence, fmt.Errorf("fixed-sigma diagnostic requires the existing Fast Q-prefix capacity checker")
	}
	evidence.Capacity, err = capacityChecker.CheckQPrefixAdditionCapacity(ringQ, ct, &candidateCT, coefficients, rows)
	if err != nil {
		return evidence, fmt.Errorf("fixed-sigma diagnostic candidate exceeds or cannot prove Q-prefix capacity: %w", err)
	}
	// Commit only c0 after all rows and centered capacity have passed. Level,
	// scale, metadata, c1, and dormant logical limbs remain untouched.
	ct.Value[0] = *candidate
	evidence.DispatchEvidence = "explicit Fast ckks.Evaluator helper; q0..q3 coefficient residues transformed per NTT subring; existing fastcore centered-Q-prefix addition-capacity check; c1/dormant limbs untouched"
	return evidence, nil
}

func signedCoefficientResidue(value int64, modulus uint64) uint64 {
	if value >= 0 {
		return uint64(value) % modulus
	}
	magnitude := uint64(-(value + 1)) + 1
	remainder := magnitude % modulus
	if remainder == 0 {
		return 0
	}
	return modulus - remainder
}
