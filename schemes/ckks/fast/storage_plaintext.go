package fast

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// FastStoragePlaintext is an immutable width-3 private-F mirror of one
// encoded CKKS plaintext polynomial. Its rows are ordinary-residue NTT
// polynomials and are physically separate from LogicalQ rows.
type FastStoragePlaintext struct {
	params   ckks.Parameters
	levelQ   int
	scale    rlwe.Scale
	metadata rlwe.PlaintextMetaData
	value    ring.Poly
	maxAbs   *big.Int
	l1Norm   *big.Int
	basis    fastStorageBasis
}

// NewFastStoragePlaintextMirror converts one complete NTT/Montgomery
// LogicalQ polynomial into its centered integer polynomial and then into the
// fixed width-3 private-F NTT representation. The input is never modified.
func NewFastStoragePlaintextMirror(params ckks.Parameters, logicalQ ring.Poly, levelQ int, metadata rlwe.PlaintextMetaData, isNTT, isMontgomery bool) (*FastStoragePlaintext, error) {
	if params.N() == 0 || params.RingType() != ring.Standard {
		return nil, errors.New("private-F plaintext mirror requires initialized Standard CKKS parameters")
	}
	if levelQ < 0 || levelQ > params.MaxLevel() {
		return nil, fmt.Errorf("plaintext LevelQ %d is outside [0,%d]", levelQ, params.MaxLevel())
	}
	if !isNTT || !isMontgomery {
		return nil, errors.New("LogicalQ plaintext mirror requires NTT and Montgomery source rows")
	}
	if metadata.Scale.Cmp(rlwe.NewScale(0)) != 1 {
		return nil, errors.New("plaintext Scale must be positive")
	}
	if len(logicalQ.Coeffs) != levelQ+1 {
		return nil, fmt.Errorf("LogicalQ plaintext has %d rows, expected %d through LevelQ %d", len(logicalQ.Coeffs), levelQ+1, levelQ)
	}
	for row, coeffs := range logicalQ.Coeffs {
		if len(coeffs) != params.N() {
			return nil, fmt.Errorf("LogicalQ plaintext row %d has length %d, expected N=%d", row, len(coeffs), params.N())
		}
		q := params.Q()[row]
		for coefficient, residue := range coeffs {
			if residue >= q {
				return nil, fmt.Errorf("LogicalQ plaintext row %d coefficient %d is not canonical modulo q%d", row, coefficient, row)
			}
		}
	}

	basis, err := fastStorageBasisForLogN(params.LogN())
	if err != nil {
		return nil, err
	}
	logicalRing := params.RingQ().AtLevel(levelQ)
	coefficients := ring.NewPoly(params.N(), levelQ)
	for row := range logicalQ.Coeffs {
		copy(coefficients.Coeffs[row], logicalQ.Coeffs[row])
	}
	logicalRing.IMForm(coefficients, coefficients)
	logicalRing.INTT(coefficients, coefficients)

	integers := make([]*big.Int, params.N())
	for i := range integers {
		integers[i] = new(big.Int)
	}
	logicalRing.PolyToBigintCentered(coefficients, 1, integers)

	maxAbs, l1Norm := new(big.Int), new(big.Int)
	fCoeff := ring.NewPoly(params.N(), 2)
	for coefficient, integer := range integers {
		magnitude := new(big.Int).Abs(new(big.Int).Set(integer))
		if magnitude.Cmp(maxAbs) > 0 {
			maxAbs.Set(magnitude)
		}
		l1Norm.Add(l1Norm, magnitude)
		unique, err := basis.hasUniqueCenteredRepresentation(integer, 3)
		if err != nil {
			return nil, err
		}
		if !unique {
			return nil, fmt.Errorf("plaintext coefficient %d magnitude %s exceeds strict width-3 centered capacity", coefficient, magnitude)
		}
		encoded, err := basis.encode(integer, 3)
		if err != nil {
			return nil, fmt.Errorf("encode plaintext coefficient %d into private F: %w", coefficient, err)
		}
		for row := 0; row < 3; row++ {
			fCoeff.Coeffs[row][coefficient] = encoded[row]
		}
	}
	for row := 0; row < 3; row++ {
		subring, _ := basis.subring(row)
		subring.NTT(fCoeff.Coeffs[row], fCoeff.Coeffs[row])
	}

	if err := verifyFastStoragePlaintextMirror(params, levelQ, logicalQ, fCoeff, basis, integers); err != nil {
		return nil, err
	}
	return &FastStoragePlaintext{
		params: params, levelQ: levelQ, scale: cloneFastScale(metadata.Scale),
		metadata: cloneFastPlaintextMetadata(metadata), value: fCoeff,
		maxAbs: new(big.Int).Set(maxAbs), l1Norm: new(big.Int).Set(l1Norm), basis: basis,
	}, nil
}

func verifyFastStoragePlaintextMirror(params ckks.Parameters, levelQ int, original, privateF ring.Poly, basis fastStorageBasis, sourceIntegers []*big.Int) error {
	fCoefficients := make([][]uint64, 3)
	for row := 0; row < 3; row++ {
		fCoefficients[row] = make([]uint64, params.N())
		subring, _ := basis.subring(row)
		subring.INTT(privateF.Coeffs[row], fCoefficients[row])
	}
	integers := make([]*big.Int, params.N())
	for coefficient := range integers {
		var residues [3]uint64
		for row := 0; row < 3; row++ {
			residues[row] = fCoefficients[row][coefficient]
		}
		integer, err := basis.decode(residues, 3)
		if err != nil {
			return fmt.Errorf("reconstruct private-F plaintext coefficient %d: %w", coefficient, err)
		}
		integers[coefficient] = integer
		if integer.Cmp(sourceIntegers[coefficient]) != 0 {
			return fmt.Errorf("private-F plaintext coefficient %d does not reconstruct the centered LogicalQ integer", coefficient)
		}
	}
	logical := ring.NewPoly(params.N(), levelQ)
	for coefficient, integer := range integers {
		for row := 0; row <= levelQ; row++ {
			q := new(big.Int).SetUint64(params.Q()[row])
			logical.Coeffs[row][coefficient] = new(big.Int).Mod(new(big.Int).Set(integer), q).Uint64()
		}
	}
	logicalRing := params.RingQ().AtLevel(levelQ)
	logicalRing.NTT(logical, logical)
	logicalRing.MForm(logical, logical)
	for row := 0; row <= levelQ; row++ {
		for coefficient, got := range logical.Coeffs[row] {
			if got != original.Coeffs[row][coefficient] {
				return fmt.Errorf("private-F plaintext mirror round trip differs from source at q%d coefficient %d", row, coefficient)
			}
		}
	}
	return nil
}

func cloneFastPlaintextMetadata(metadata rlwe.PlaintextMetaData) rlwe.PlaintextMetaData {
	copy := metadata
	copy.Scale = cloneFastScale(metadata.Scale)
	return copy
}

// CopyNew returns an independent mirror, including its exact bound metadata.
func (pt *FastStoragePlaintext) CopyNew() *FastStoragePlaintext {
	if pt == nil {
		return nil
	}
	cloned := *pt
	cloned.scale = cloneFastScale(pt.scale)
	cloned.metadata = cloneFastPlaintextMetadata(pt.metadata)
	cloned.maxAbs = new(big.Int).Set(pt.maxAbs)
	cloned.l1Norm = new(big.Int).Set(pt.l1Norm)
	cloned.value = ring.NewPoly(pt.params.N(), 2)
	for row := range pt.value.Coeffs {
		copy(cloned.value.Coeffs[row], pt.value.Coeffs[row])
	}
	return &cloned
}

func (pt *FastStoragePlaintext) LevelQ() int {
	if pt == nil {
		return -1
	}
	return pt.levelQ
}

func (pt *FastStoragePlaintext) Scale() rlwe.Scale {
	if pt == nil {
		return rlwe.Scale{}
	}
	return cloneFastScale(pt.scale)
}

func (pt *FastStoragePlaintext) MetaData() rlwe.PlaintextMetaData {
	if pt == nil {
		return rlwe.PlaintextMetaData{}
	}
	return cloneFastPlaintextMetadata(pt.metadata)
}

func (pt *FastStoragePlaintext) MaxAbsCoefficient() *big.Int {
	if pt == nil || pt.maxAbs == nil {
		return nil
	}
	return new(big.Int).Set(pt.maxAbs)
}

func (pt *FastStoragePlaintext) L1Norm() *big.Int {
	if pt == nil || pt.l1Norm == nil {
		return nil
	}
	return new(big.Int).Set(pt.l1Norm)
}

func (pt *FastStoragePlaintext) StorageWidth() int {
	if pt == nil {
		return 0
	}
	return 3
}

func (pt *FastStoragePlaintext) IsNTT() bool { return pt != nil }

func (pt *FastStoragePlaintext) IsMontgomery() bool { return false }

func validateFastStoragePlaintext(pt *FastStoragePlaintext) error {
	if pt == nil {
		return errors.New("Fast storage plaintext cannot be nil")
	}
	if pt.params.N() == 0 || pt.params.RingType() != ring.Standard || pt.levelQ < 0 || pt.levelQ > pt.params.MaxLevel() {
		return errors.New("Fast storage plaintext has invalid CKKS parameters or LevelQ")
	}
	if pt.basis.logN != pt.params.LogN() {
		return errors.New("Fast storage plaintext basis does not match CKKS ring degree")
	}
	if pt.scale.Cmp(rlwe.NewScale(0)) != 1 || pt.maxAbs == nil || pt.l1Norm == nil || pt.maxAbs.Sign() < 0 || pt.l1Norm.Sign() < 0 {
		return errors.New("Fast storage plaintext has invalid Scale or coefficient bounds")
	}
	if len(pt.value.Coeffs) != 3 {
		return fmt.Errorf("Fast storage plaintext has %d F rows, expected width 3", len(pt.value.Coeffs))
	}
	for row := 0; row < 3; row++ {
		subring, err := pt.basis.subring(row)
		if err != nil {
			return err
		}
		if len(pt.value.Coeffs[row]) != pt.params.N() {
			return fmt.Errorf("Fast storage plaintext F%d row length mismatch", row)
		}
		for coefficient, residue := range pt.value.Coeffs[row] {
			if residue >= subring.Modulus {
				return fmt.Errorf("Fast storage plaintext F%d coefficient %d is not canonical", row, coefficient)
			}
		}
	}
	return nil
}
