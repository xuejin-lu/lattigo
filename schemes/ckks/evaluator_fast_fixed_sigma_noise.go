package ckks

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/internal/fastcore"
)

// FixedSigmaCoefficientNoiseEvidence describes one explicit diagnostic-only
// fixed-sigma injection. Coefficients are rounded signed Gaussian integers;
// the same signed polynomial is reduced into every maintained Q row.
type FixedSigmaCoefficientNoiseEvidence struct {
	Sigma               float64                           `json:"sigma_coefficient"`
	Seed                int64                             `json:"sampling_seed"`
	Level               int                               `json:"level"`
	Scale               string                            `json:"scale"`
	NoOp                bool                              `json:"no_op"`
	SampleCount         int                               `json:"coefficient_sample_count"`
	CoefficientMean     float64                           `json:"coefficient_mean"`
	CoefficientVariance float64                           `json:"coefficient_variance"`
	CoefficientRMS      float64                           `json:"coefficient_rms"`
	CoefficientMaxAbs   uint64                            `json:"coefficient_max_abs"`
	CoefficientSHA256   string                            `json:"coefficient_sha256"`
	SamplingAlgorithm   string                            `json:"sampling_algorithm"`
	Capacity            fastcore.QPrefixCapacityEvidence  `json:"qprefix_capacity"`
	DispatchEvidence    string                            `json:"dispatch_evidence"`
	Phases              []FixedSigmaCoefficientNoisePhase `json:"phase_measurements,omitempty"`
}

// FixedSigmaCoefficientNoisePhase records a bounded helper subphase for CP-D.
type FixedSigmaCoefficientNoisePhase struct {
	Name                 string `json:"name"`
	ElapsedNS            int64  `json:"elapsed_ns"`
	AllocBytes           uint64 `json:"alloc_bytes"`
	Allocs               uint64 `json:"allocs"`
	MemoryStatsCollected bool   `json:"memory_stats_collected"`
}

// InjectFixedSigmaCoefficientNoiseDiagnostic explicitly adds one sampled
// coefficient-domain polynomial to a Fast zero-secret ciphertext's c0. It is
// diagnostic-only, supports the isolated LogN13 task boundary at Level 4 or
// Level 3 with a four-row Q-prefix, and is never called by ordinary CKKS
// operations.
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
	if params.LogN() != 13 || params.N() != 1<<13 || (ct.Level() != 4 && ct.Level() != 3) || !ct.Scale.Equal(rlwe.NewScale(math.Exp2(45))) ||
		!ct.IsNTT || ct.IsMontgomery || !eval.GetRLWEParameters().NTTFlag() {
		return evidence, fmt.Errorf("fixed-sigma diagnostic requires the LogN13, Level4 or Level3, Scale2^45 NTT non-Montgomery output profile")
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

	return eval.injectFixedSigmaCoefficientNoise(ct, sigma, seed, rows, true)
}

// InjectFixedSigmaBootstrapMulNoiseDiagnostic injects the same coefficient-
// domain diagnostic noise as InjectFixedSigmaCoefficientNoiseDiagnostic, but
// accepts every active Q-prefix level and positive Scale used by the bounded
// LogN13 Bootstrap multiplication path. It remains explicit and is never
// called by ordinary CKKS operations. The historical post-Mul/Rescale helper
// above retains its Level3/4 and Scale2^45 contract.
func (eval Evaluator) InjectFixedSigmaBootstrapMulNoiseDiagnostic(ct *rlwe.Ciphertext, sigma float64, seed int64) (FixedSigmaCoefficientNoiseEvidence, error) {
	evidence := FixedSigmaCoefficientNoiseEvidence{Sigma: sigma, Seed: seed}
	if math.IsNaN(sigma) || math.IsInf(sigma, 0) || sigma <= 0 {
		return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul diagnostic requires finite sigma > 0")
	}
	if eval.Encoder == nil || eval.Evaluator == nil {
		return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul diagnostic requires an initialized CKKS evaluator")
	}
	if !eval.usesFastCKKSZeroSecretSimulation() {
		return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul diagnostic is available only for the Fast zero-secret CKKS capability")
	}
	params := eval.GetParameters()
	if params.LogN() != 13 || params.N() != 1<<13 {
		return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul diagnostic requires the frozen LogN13 profile")
	}
	if _, err := validateFastCKKSZeroSecretMulShape("Bootstrap Mul target", ct, params); err != nil {
		return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul target: %w", err)
	}
	if ct.Scale.Value.Sign() <= 0 || !ct.IsNTT || ct.IsMontgomery || !eval.GetRLWEParameters().NTTFlag() {
		return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul target requires a positive Scale and NTT non-Montgomery representation")
	}
	rows, err := fastcore.QPrefixWidth(ct.Level())
	if err != nil {
		return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul Q-prefix width: %w", err)
	}
	ringQ := params.RingQ()
	if err := fastcore.ValidatePrefixRows(ringQ, ct.Level(), rows, ct.Value...); err != nil {
		return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul Q-prefix backing: %w", err)
	}
	for row := 0; row < rows; row++ {
		modulus := ringQ.SubRings[row].Modulus
		for coefficient, value := range ct.Value[0].Coeffs[row] {
			if value >= modulus {
				return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul c0 q%d coefficient %d is not a canonical residue", row, coefficient)
			}
		}
		for coefficient, value := range ct.Value[1].Coeffs[row] {
			if value != 0 {
				return evidence, fmt.Errorf("fixed-sigma Bootstrap Mul diagnostic requires c1=0; q%d coefficient %d is nonzero", row, coefficient)
			}
		}
	}
	return eval.injectFixedSigmaCoefficientNoise(ct, sigma, seed, rows, false)
}

func (eval Evaluator) injectFixedSigmaCoefficientNoise(ct *rlwe.Ciphertext, sigma float64, seed int64, rows int, collectMemoryStats bool) (FixedSigmaCoefficientNoiseEvidence, error) {
	params := eval.GetParameters()
	ringQ := params.RingQ()
	evidence := FixedSigmaCoefficientNoiseEvidence{Sigma: sigma, Seed: seed, Level: ct.Level(), Scale: ct.Scale.Value.Text('g', -1)}
	n := params.N()
	evidence.SamplingAlgorithm = "Go math/rand.NewSource(seed).NormFloat64(), rounded to nearest integer with math.Round (ties away from zero)"
	var coefficients []int64
	phase, phaseErr := measureFixedSigmaCoefficientNoisePhase("sampling_hash_moments", collectMemoryStats, func() error {
		coefficients = make([]int64, n)
		rng := rand.New(rand.NewSource(seed))
		hasher := sha256.New()
		var mean, m2, sumSquares float64
		var maxAbs uint64
		var encoded [8]byte
		for index := range coefficients {
			sampled := math.Round(rng.NormFloat64() * sigma)
			if math.IsNaN(sampled) || math.IsInf(sampled, 0) || sampled >= math.Exp2(63) || sampled < -math.Exp2(63) {
				return fmt.Errorf("fixed-sigma diagnostic sample %d cannot be represented as int64", index)
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
		evidence.SampleCount = n
		evidence.CoefficientMean = mean
		evidence.CoefficientVariance = m2 / float64(n)
		evidence.CoefficientRMS = math.Sqrt(sumSquares / float64(n))
		evidence.CoefficientMaxAbs = maxAbs
		evidence.CoefficientSHA256 = hex.EncodeToString(hasher.Sum(nil))
		return nil
	})
	evidence.Phases = append(evidence.Phases, phase)
	if phaseErr != nil {
		return evidence, phaseErr
	}

	var noiseNTT ring.Poly
	phase, phaseErr = measureFixedSigmaCoefficientNoisePhase("coefficient_to_qprefix_ntt", collectMemoryStats, func() error {
		noiseNTT = ring.NewPoly(n, rows-1)
		for row := 0; row < rows; row++ {
			modulus := ringQ.SubRings[row].Modulus
			for coefficient, value := range coefficients {
				noiseNTT.Coeffs[row][coefficient] = signedCoefficientResidue(value, modulus)
			}
			ringQ.SubRings[row].NTT(noiseNTT.Coeffs[row], noiseNTT.Coeffs[row])
		}
		return nil
	})
	evidence.Phases = append(evidence.Phases, phase)
	if phaseErr != nil {
		return evidence, phaseErr
	}

	var candidate *ring.Poly
	var candidateCT rlwe.Ciphertext
	phase, phaseErr = measureFixedSigmaCoefficientNoisePhase("candidate_c0_add", collectMemoryStats, func() error {
		candidate = ct.Value[0].CopyNew()
		for row := 0; row < rows; row++ {
			ringQ.SubRings[row].Add(candidate.Coeffs[row], noiseNTT.Coeffs[row], candidate.Coeffs[row])
		}
		candidateCT = *ct
		candidateCT.Value = append([]ring.Poly(nil), ct.Value...)
		candidateCT.Value[0] = *candidate
		return nil
	})
	evidence.Phases = append(evidence.Phases, phase)
	if phaseErr != nil {
		return evidence, phaseErr
	}
	capacityChecker, ok := eval.rescaleCore.(interface {
		CheckQPrefixAdditionCapacity(*ring.Ring, *rlwe.Ciphertext, *rlwe.Ciphertext, []int64, int) (fastcore.QPrefixCapacityEvidence, error)
	})
	if !ok {
		return evidence, fmt.Errorf("fixed-sigma diagnostic requires the existing Fast Q-prefix capacity checker")
	}
	phase, phaseErr = measureFixedSigmaCoefficientNoisePhase("qprefix_centered_capacity_check", collectMemoryStats, func() error {
		var checkErr error
		evidence.Capacity, checkErr = capacityChecker.CheckQPrefixAdditionCapacity(ringQ, ct, &candidateCT, coefficients, rows)
		return checkErr
	})
	evidence.Phases = append(evidence.Phases, phase)
	if phaseErr != nil {
		return evidence, fmt.Errorf("fixed-sigma diagnostic candidate exceeds or cannot prove Q-prefix capacity: %w", phaseErr)
	}
	// Commit only c0 after all rows and centered capacity have passed. Level,
	// scale, metadata, c1, and dormant logical limbs remain untouched.
	phase, phaseErr = measureFixedSigmaCoefficientNoisePhase("commit_c0", collectMemoryStats, func() error {
		ct.Value[0] = *candidate
		return nil
	})
	evidence.Phases = append(evidence.Phases, phase)
	if phaseErr != nil {
		return evidence, phaseErr
	}
	evidence.DispatchEvidence = "explicit Fast ckks.Evaluator helper; q0..q3 coefficient residues transformed per NTT subring; existing fastcore centered-Q-prefix addition-capacity check; c1/dormant limbs untouched"
	return evidence, nil
}

func measureFixedSigmaCoefficientNoisePhase(name string, collectMemoryStats bool, operation func() error) (FixedSigmaCoefficientNoisePhase, error) {
	var before, after runtime.MemStats
	if collectMemoryStats {
		runtime.ReadMemStats(&before)
	}
	started := time.Now()
	err := operation()
	elapsed := time.Since(started).Nanoseconds()
	phase := FixedSigmaCoefficientNoisePhase{Name: name, ElapsedNS: elapsed, MemoryStatsCollected: collectMemoryStats}
	if collectMemoryStats {
		runtime.ReadMemStats(&after)
		phase.AllocBytes = after.TotalAlloc - before.TotalAlloc
		phase.Allocs = after.Mallocs - before.Mallocs
	}
	return phase, err
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
