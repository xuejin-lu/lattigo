package fast

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// FastCiphertext keeps private-F residue polynomials physically separate from
// ordinary LogicalQ ciphertext rows. Logical level and storage width are
// independent: each Value polynomial has exactly activeStorageWidth rows.
// The fields are private so callers cannot relabel storage rows or mutate the
// representation flags without crossing an explicit conversion boundary.
type FastCiphertext struct {
	params             ckks.Parameters
	logicalLevel       int
	activeStorageWidth int
	value              []ring.Poly
	metadata           rlwe.MetaData
	basis              fastStorageBasis
}

// NewFastCiphertext allocates a zero Fast ciphertext with an explicit logical
// level and independent private storage width.
func NewFastCiphertext(params ckks.Parameters, degree, logicalLevel, storageWidth int) (*FastCiphertext, error) {
	if err := validateFastCiphertextParameters(params, degree, logicalLevel, storageWidth); err != nil {
		return nil, err
	}
	basis, err := newFastStorageBasis(params.LogN())
	if err != nil {
		return nil, err
	}
	return newFastCiphertextWithBasis(params, degree, logicalLevel, storageWidth, basis), nil
}

func newFastCiphertextWithBasis(params ckks.Parameters, degree, logicalLevel, storageWidth int, basis fastStorageBasis) *FastCiphertext {
	n := params.N()
	value := make([]ring.Poly, degree+1)
	for i := range value {
		value[i] = ring.NewPoly(n, storageWidth-1)
	}
	return &FastCiphertext{
		params:             params,
		logicalLevel:       logicalLevel,
		activeStorageWidth: storageWidth,
		value:              value,
		metadata: rlwe.MetaData{
			PlaintextMetaData: rlwe.PlaintextMetaData{
				Scale:         cloneFastScale(params.DefaultScale()),
				LogDimensions: params.LogMaxDimensions(),
				IsBatched:     true,
			},
			CiphertextMetaData: rlwe.CiphertextMetaData{
				IsNTT: params.GetRLWEParameters().NTTFlag(),
			},
		},
		basis: basis,
	}
}

func validateFastCiphertextParameters(params ckks.Parameters, degree, logicalLevel, storageWidth int) error {
	if params.N() == 0 {
		return fmt.Errorf("CKKS parameters cannot be empty")
	}
	if params.RingType() != ring.Standard {
		return fmt.Errorf("Fast storage ciphertext requires Standard ring parameters")
	}
	if degree < 0 {
		return fmt.Errorf("ciphertext degree must be non-negative: %d", degree)
	}
	if logicalLevel < 0 || logicalLevel > params.MaxLevel() {
		return fmt.Errorf("logical level %d is outside [0,%d]", logicalLevel, params.MaxLevel())
	}
	return validateStorageWidth(storageWidth)
}

// CopyNew returns a deep copy of the value rows and metadata. Immutable CKKS
// parameters and NTT subring tables are shared by value/reference.
func (ct *FastCiphertext) CopyNew() *FastCiphertext {
	if ct == nil {
		return nil
	}
	copy := *ct
	copy.metadata = cloneFastMetadata(ct.metadata)
	copy.value = make([]ring.Poly, len(ct.value))
	for i := range ct.value {
		copy.value[i].Coeffs = make([][]uint64, len(ct.value[i].Coeffs))
		for j := range ct.value[i].Coeffs {
			copy.value[i].Coeffs[j] = append([]uint64(nil), ct.value[i].Coeffs[j]...)
		}
	}
	return &copy
}

// ResizeDegree changes only the number of ciphertext components. Existing
// rows are retained; new components are zero-filled in the current F basis.
func (ct *FastCiphertext) ResizeDegree(degree int) error {
	if ct == nil {
		return fmt.Errorf("Fast ciphertext cannot be nil")
	}
	if degree < 0 {
		return fmt.Errorf("ciphertext degree must be non-negative: %d", degree)
	}
	if ct.N() == 0 {
		return fmt.Errorf("Fast ciphertext is not initialized")
	}
	if err := validateStorageWidth(ct.activeStorageWidth); err != nil {
		return err
	}
	wantComponents := degree + 1
	if wantComponents < len(ct.value) {
		ct.value = ct.value[:wantComponents]
		return nil
	}
	for len(ct.value) < wantComponents {
		poly := ring.NewPoly(ct.N(), ct.activeStorageWidth-1)
		ct.value = append(ct.value, poly)
	}
	return nil
}

// SetLogicalLevel changes CKKS level metadata only. It never resizes or
// rewrites the private-F residue rows.
func (ct *FastCiphertext) SetLogicalLevel(level int) error {
	if ct == nil {
		return fmt.Errorf("Fast ciphertext cannot be nil")
	}
	if level < 0 || level > ct.params.MaxLevel() {
		return fmt.Errorf("logical level %d is outside [0,%d]", level, ct.params.MaxLevel())
	}
	ct.logicalLevel = level
	return nil
}

// SetScale updates the CKKS scale metadata without changing the storage rows.
func (ct *FastCiphertext) SetScale(scale rlwe.Scale) error {
	if ct == nil {
		return fmt.Errorf("Fast ciphertext cannot be nil")
	}
	if scale.Cmp(rlwe.NewScale(0)) != 1 {
		return fmt.Errorf("Fast ciphertext scale must be positive")
	}
	ct.metadata.Scale = cloneFastScale(scale)
	return nil
}

func (ct *FastCiphertext) LogicalLevel() int {
	if ct == nil {
		return -1
	}
	return ct.logicalLevel
}

func (ct *FastCiphertext) StorageWidth() int {
	if ct == nil {
		return 0
	}
	return ct.activeStorageWidth
}

func (ct *FastCiphertext) N() int {
	if ct == nil {
		return 0
	}
	return ct.params.N()
}

func (ct *FastCiphertext) LogN() int {
	if ct == nil || ct.params.N() == 0 {
		return 0
	}
	return ct.params.LogN()
}

func (ct *FastCiphertext) Degree() int {
	if ct == nil {
		return -1
	}
	return len(ct.value) - 1
}

func (ct *FastCiphertext) Scale() rlwe.Scale {
	if ct == nil {
		return rlwe.Scale{}
	}
	return cloneFastScale(ct.metadata.Scale)
}

func (ct *FastCiphertext) IsNTT() bool {
	return ct != nil && ct.metadata.IsNTT
}

func (ct *FastCiphertext) IsMontgomery() bool {
	return ct != nil && ct.metadata.IsMontgomery
}

func (ct *FastCiphertext) MetaData() rlwe.MetaData {
	if ct == nil {
		return rlwe.MetaData{}
	}
	return cloneFastMetadata(ct.metadata)
}

func (ct *FastCiphertext) Parameters() ckks.Parameters {
	if ct == nil {
		return ckks.Parameters{}
	}
	return ct.params
}

func cloneFastScale(scale rlwe.Scale) rlwe.Scale {
	copy := rlwe.Scale{Value: *new(big.Float).Copy(&scale.Value)}
	if scale.Mod != nil {
		copy.Mod = new(big.Int).Set(scale.Mod)
	}
	return copy
}

func cloneFastMetadata(metadata rlwe.MetaData) rlwe.MetaData {
	copy := metadata
	copy.Scale = cloneFastScale(metadata.Scale)
	return copy
}
