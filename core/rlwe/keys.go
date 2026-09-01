package rlwe

import (
	"bufio"
	"fmt"
	"io"
	"slices"

	"github.com/google/go-cmp/cmp"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils/buffer"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
	"github.com/tuneinsight/lattigo/v6/utils/structs"
)

// KeyLayout explicitly defines whether a key incorporates the deterministic
// key-dependent term (Standard) or is purely an error sample (Fast).
type KeyLayout int

const (
	KeyLayoutStandard KeyLayout = 0
	KeyLayoutFast     KeyLayout = 1
)

type SecretKey struct {
	Value ringqp.Poly
}

func NewSecretKey(params ParameterProvider) *SecretKey {
	return &SecretKey{Value: params.GetRLWEParameters().RingQP().NewPoly()}
}

func (sk SecretKey) Equal(other *SecretKey) bool {
	return cmp.Equal(sk.Value, other.Value)
}

func (sk SecretKey) LevelQ() int {
	return sk.Value.Q.Level()
}

func (sk SecretKey) LevelP() int {
	return sk.Value.P.Level()
}

func (sk SecretKey) CopyNew() *SecretKey {
	return &SecretKey{*sk.Value.CopyNew()}
}

func (sk SecretKey) BinarySize() (dataLen int) {
	return sk.Value.BinarySize()
}

func (sk SecretKey) WriteTo(w io.Writer) (n int64, err error) {
	return sk.Value.WriteTo(w)
}

func (sk *SecretKey) ReadFrom(r io.Reader) (n int64, err error) {
	return sk.Value.ReadFrom(r)
}

func (sk SecretKey) MarshalBinary() (p []byte, err error) {
	return sk.Value.MarshalBinary()
}

func (sk *SecretKey) UnmarshalBinary(p []byte) (err error) {
	return sk.Value.UnmarshalBinary(p)
}

func (sk *SecretKey) isEncryptionKey() {}

type VectorQP []ringqp.Poly

func NewVectorQP(params ParameterProvider, size, levelQ, levelP int) (v VectorQP) {
	rqp := params.GetRLWEParameters().RingQP().AtLevel(levelQ, levelP)
	v = make(VectorQP, size)
	for i := range v {
		v[i] = rqp.NewPoly()
	}
	return
}

func (p VectorQP) LevelQ() int {
	if len(p) == 0 {
		return -1
	}
	return p[0].LevelQ()
}

func (p VectorQP) LevelP() int {
	if len(p) == 0 {
		return -1
	}
	return p[0].LevelP()
}

func (p VectorQP) CopyNew() *VectorQP {
	m := make([]ringqp.Poly, len(p))
	for i := range p {
		m[i] = *p[i].CopyNew()
	}
	v := VectorQP(m)
	return &v
}

func (p VectorQP) Equal(other *VectorQP) (equal bool) {
	if len(p) != len(*other) {
		return false
	}
	equal = true
	for i := range p {
		equal = equal && p[i].Equal(&(*other)[i])
	}
	return
}

func (p VectorQP) BinarySize() int {
	return structs.Vector[ringqp.Poly](p[:]).BinarySize()
}

func (p VectorQP) WriteTo(w io.Writer) (n int64, err error) {
	v := structs.Vector[ringqp.Poly](p[:])
	return v.WriteTo(w)
}

func (p *VectorQP) ReadFrom(r io.Reader) (n int64, err error) {
	v := structs.Vector[ringqp.Poly](*p)
	n, err = v.ReadFrom(r)
	*p = VectorQP(v)
	return
}

func (p VectorQP) MarshalBinary() ([]byte, error) {
	buf := buffer.NewBufferSize(p.BinarySize())
	_, err := p.WriteTo(buf)
	return buf.Bytes(), err
}

func (p *VectorQP) UnmarshalBinary(b []byte) error {
	_, err := p.ReadFrom(buffer.NewBuffer(b))
	return err
}

type PublicKey struct {
	Layout KeyLayout
	Value  VectorQP
}

func NewPublicKey(params ParameterProvider) (pk *PublicKey) {
	p := params.GetRLWEParameters()
	return &PublicKey{Layout: KeyLayoutStandard, Value: NewVectorQP(params, 2, p.MaxLevelQ(), p.MaxLevelP())}
}

func (p PublicKey) LevelQ() int {
	return p.Value.LevelQ()
}

func (p PublicKey) LevelP() int {
	return p.Value.LevelP()
}

func (p PublicKey) CopyNew() *PublicKey {
	return &PublicKey{Layout: p.Layout, Value: *p.Value.CopyNew()}
}

func (p PublicKey) Equal(other *PublicKey) bool {
	return p.Value.Equal(&other.Value)
}

func (p PublicKey) BinarySize() int {
	return p.Value.BinarySize()
}

func (p PublicKey) WriteTo(w io.Writer) (n int64, err error) {
	if p.Layout == KeyLayoutFast {
		return 0, fmt.Errorf("unsupported: Fast keys cannot be serialized")
	}
	return p.Value.WriteTo(w)
}

func (p *PublicKey) ReadFrom(r io.Reader) (n int64, err error) {
	return p.Value.ReadFrom(r)
}

func (p PublicKey) MarshalBinary() ([]byte, error) {
	if p.Layout == KeyLayoutFast {
		return nil, fmt.Errorf("unsupported: Fast keys cannot be serialized")
	}
	return p.Value.MarshalBinary()
}

func (p *PublicKey) UnmarshalBinary(b []byte) error {
	return p.Value.UnmarshalBinary(b)
}

func (p *PublicKey) isEncryptionKey() {}

type EvaluationKey struct {
	GadgetCiphertext
	Seed   *[32]byte
	Layout KeyLayout
}

type EvaluationKeyParameters struct {
	LevelQ               *int
	LevelP               *int
	BaseTwoDecomposition *int
	Compressed           bool
}

func ResolveEvaluationKeyParameters(params Parameters, evkParams []EvaluationKeyParameters) (levelQ, levelP, BaseTwoDecomposition int, compressed bool) {
	if len(evkParams) != 0 {
		if evkParams[0].LevelQ == nil {
			levelQ = params.MaxLevelQ()
		} else {
			levelQ = *evkParams[0].LevelQ
		}

		if evkParams[0].LevelP == nil {
			levelP = params.MaxLevelP()
		} else {
			levelP = *evkParams[0].LevelP
		}

		if evkParams[0].BaseTwoDecomposition != nil {
			BaseTwoDecomposition = *evkParams[0].BaseTwoDecomposition
		}
		compressed = evkParams[0].Compressed
	} else {
		levelQ = params.MaxLevelQ()
		levelP = params.MaxLevelP()
	}
	return
}

func NewEvaluationKey(params ParameterProvider, evkParams ...EvaluationKeyParameters) *EvaluationKey {
	p := *params.GetRLWEParameters()
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(p, evkParams)
	return newEvaluationKey(p, levelQ, levelP, BaseTwoDecomposition, compressed)
}

func newEvaluationKey(params Parameters, levelQ, levelP, BaseTwoDecomposition int, compressed bool) *EvaluationKey {
	degree := 1
	if compressed {
		degree = 0
	}
	return &EvaluationKey{
		Layout:           KeyLayoutStandard,
		GadgetCiphertext: *NewGadgetCiphertext(params, degree, levelQ, levelP, BaseTwoDecomposition),
	}
}

func (evk EvaluationKey) IsCompressed() bool {
	return evk.Degree() == 0
}

func (evk EvaluationKey) Expand(params ParameterProvider, buffer *GadgetCiphertext) error {
	if evk.Layout == KeyLayoutFast {
		return fmt.Errorf("unsupported: Fast EvaluationKeys cannot be expanded (explicit error material cannot be seeded)")
	}

	if !evk.IsCompressed() {
		return fmt.Errorf("evaluation key is not compressed")
	}

	if evk.Seed == nil {
		return fmt.Errorf("seed is missing")
	}

	prng, err := sampling.NewKeyedPRNG((*evk.Seed)[:])
	if err != nil {
		panic(fmt.Errorf("sampling.NewKeyedPRNG: %s", err))
	}

	levelQ := evk.LevelQ()
	levelP := evk.LevelP()

	uniformRingQPSampler := ringqp.NewUniformSampler(prng, *params.GetRLWEParameters().RingQP()).AtLevel(levelQ, levelP)
	BaseRNSDecompositionVectorSize := evk.BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := evk.BaseTwoDecompositionVectorSize()

	if buffer != nil {
		if have := buffer.Degree(); have != 0 {
			return fmt.Errorf("invalid buffer, degree should be 0 but is %d", have)
		}
		if have := buffer.BaseRNSDecompositionVectorSize(); have != BaseRNSDecompositionVectorSize {
			return fmt.Errorf("invalid buffer BaseRNSDecompositionVectorSize, should be %d but is %d", have, BaseRNSDecompositionVectorSize)
		}
		if have := buffer.BaseTwoDecompositionVectorSize(); !slices.Equal(have, BaseTwoDecompositionVectorSize) {
			return fmt.Errorf("invalid buffer BaseTwoDecompositionVectorSize, should be %v but is %v", have, BaseTwoDecompositionVectorSize)
		}
		if have := buffer.LevelQ(); have != levelQ {
			return fmt.Errorf("invalid buffer levelQ, should be %d but is %d", levelQ, have)
		}
		if have := buffer.LevelP(); have != levelP {
			return fmt.Errorf("invalid buffer levelP, should be %d but is %d", levelP, have)
		}

	} else {
		buffer = NewGadgetCiphertext(params, 0, levelQ, levelP, evk.BaseTwoDecomposition)
	}

	value := make(structs.Matrix[VectorQP], BaseRNSDecompositionVectorSize)
	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		value[i] = make([]VectorQP, BaseTwoDecompositionVectorSize[i])
		for j := range value[i] {
			uniformRingQPSampler.Read(buffer.Value[i][j][0])
			evk.Value[i][j] = VectorQP{evk.Value[i][j][0], buffer.Value[i][j][0]}
		}
	}

	return nil
}

func (evk EvaluationKey) BinarySize() (size int) {
	if evk.Seed != nil {
		return evk.GadgetCiphertext.BinarySize() + len(*evk.Seed)
	}
	return evk.GadgetCiphertext.BinarySize()
}

func (evk EvaluationKey) WriteTo(w io.Writer) (n int64, err error) {
	if evk.Layout == KeyLayoutFast {
		return 0, fmt.Errorf("unsupported: Fast keys cannot be serialized")
	}

	switch w := w.(type) {
	case buffer.Writer:
		var inc int64
		if inc, err = evk.GadgetCiphertext.WriteTo(w); err != nil {
			return n + inc, err
		}
		n += inc

		if evk.IsCompressed() {
			if evk.Seed == nil {
				return n + inc, fmt.Errorf("writing compressed evaluation key: the seed is nil")
			}
			if inc, err = buffer.Write(w, (*evk.Seed)[:]); err != nil {
				return n + inc, err
			}
			n += inc
		}

		if err = w.Flush(); err != nil {
			return n, err
		}
		return

	default:
		return evk.WriteTo(bufio.NewWriter(w))
	}
}

func (evk *EvaluationKey) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:
		var inc int64
		if inc, err = evk.GadgetCiphertext.ReadFrom(r); err != nil {
			return n + inc, err
		}
		n += inc

		if evk.IsCompressed() {
			var seed [32]byte
			if inc, err = buffer.Read(r, seed[:]); err != nil {
				return n + inc, err
			}
			evk.Seed = &seed
			n += inc
		}
		return
	default:
		return evk.ReadFrom(bufio.NewReader(r))
	}
}

func (evk EvaluationKey) MarshalBinary() (p []byte, err error) {
	if evk.Layout == KeyLayoutFast {
		return nil, fmt.Errorf("unsupported: Fast keys cannot be serialized")
	}
	buf := buffer.NewBufferSize(evk.BinarySize())
	_, err = evk.WriteTo(buf)
	return buf.Bytes(), err
}

func (evk *EvaluationKey) UnmarshalBinary(p []byte) (err error) {
	_, err = evk.ReadFrom(buffer.NewBuffer(p))
	return
}

func (evk EvaluationKey) CopyNew() *EvaluationKey {
	return &EvaluationKey{Layout: evk.Layout, GadgetCiphertext: *evk.GadgetCiphertext.CopyNew()}
}

type RelinearizationKey struct {
	EvaluationKey
}

func NewRelinearizationKey(params ParameterProvider, evkParams ...EvaluationKeyParameters) *RelinearizationKey {
	p := *params.GetRLWEParameters()
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(p, evkParams)
	return newRelinearizationKey(p, levelQ, levelP, BaseTwoDecomposition, compressed)
}

func newRelinearizationKey(params Parameters, levelQ, levelP, BaseTwoDecomposition int, compressed bool) *RelinearizationKey {
	degree := 1
	if compressed {
		degree = 0
	}
	return &RelinearizationKey{EvaluationKey: EvaluationKey{Layout: KeyLayoutStandard, GadgetCiphertext: *NewGadgetCiphertext(params, degree, levelQ, levelP, BaseTwoDecomposition)}}
}

func (rlk RelinearizationKey) CopyNew() *RelinearizationKey {
	return &RelinearizationKey{EvaluationKey: *rlk.EvaluationKey.CopyNew()}
}

type GaloisKey struct {
	GaloisElement uint64
	NthRoot       uint64
	EvaluationKey
}

func NewGaloisKey(params ParameterProvider, evkParams ...EvaluationKeyParameters) *GaloisKey {
	p := *params.GetRLWEParameters()
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(p, evkParams)
	return newGaloisKey(p, levelQ, levelP, BaseTwoDecomposition, compressed)
}

func newGaloisKey(params Parameters, levelQ, levelP, BaseTwoDecomposition int, compressed bool) *GaloisKey {
	degree := 1
	if compressed {
		degree = 0
	}
	return &GaloisKey{
		EvaluationKey: EvaluationKey{
			Layout:           KeyLayoutStandard,
			GadgetCiphertext: *NewGadgetCiphertext(params, degree, levelQ, levelP, BaseTwoDecomposition),
		},
		NthRoot: params.GetRLWEParameters().RingQ().NthRoot(),
	}
}

func (gk GaloisKey) CopyNew() *GaloisKey {
	return &GaloisKey{
		GaloisElement: gk.GaloisElement,
		NthRoot:       gk.NthRoot,
		EvaluationKey: *gk.EvaluationKey.CopyNew(),
	}
}

func (gk GaloisKey) BinarySize() (size int) {
	return gk.EvaluationKey.BinarySize() + 16
}

func (gk GaloisKey) WriteTo(w io.Writer) (n int64, err error) {
	if gk.Layout == KeyLayoutFast {
		return 0, fmt.Errorf("unsupported: Fast keys cannot be serialized")
	}

	switch w := w.(type) {
	case buffer.Writer:
		var inc int64
		if inc, err = buffer.WriteUint64(w, gk.GaloisElement); err != nil {
			return n + inc, err
		}
		n += inc

		if inc, err = buffer.WriteUint64(w, gk.NthRoot); err != nil {
			return n + inc, err
		}
		n += inc

		if inc, err = gk.EvaluationKey.WriteTo(w); err != nil {
			return n + inc, err
		}
		n += inc
		return
	default:
		return gk.WriteTo(bufio.NewWriter(w))
	}
}

func (gk *GaloisKey) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:
		var inc int64
		if inc, err = buffer.ReadUint64(r, &gk.GaloisElement); err != nil {
			return n + inc, err
		}
		n += inc

		if inc, err = buffer.ReadUint64(r, &gk.NthRoot); err != nil {
			return n + inc, err
		}
		n += inc

		if inc, err = gk.EvaluationKey.ReadFrom(r); err != nil {
			return n + inc, err
		}
		n += inc
		return
	default:
		return gk.ReadFrom(bufio.NewReader(r))
	}
}

func (gk GaloisKey) MarshalBinary() (p []byte, err error) {
	if gk.Layout == KeyLayoutFast {
		return nil, fmt.Errorf("unsupported: Fast keys cannot be serialized")
	}
	buf := buffer.NewBufferSize(gk.BinarySize())
	_, err = gk.WriteTo(buf)
	return buf.Bytes(), err
}

func (gk *GaloisKey) UnmarshalBinary(p []byte) (err error) {
	_, err = gk.ReadFrom(buffer.NewBuffer(p))
	return
}

type EvaluationKeySet interface {
	GetGaloisKey(galEl uint64) (evk *GaloisKey, err error)
	GetGaloisKeysList() (galEls []uint64)
	GetRelinearizationKey() (evk *RelinearizationKey, err error)
}

type MemEvaluationKeySet struct {
	RelinearizationKey *RelinearizationKey
	GaloisKeys         structs.Map[uint64, GaloisKey]
}

func NewMemEvaluationKeySet(relinKey *RelinearizationKey, galoisKeys ...*GaloisKey) (eks *MemEvaluationKeySet) {
	eks = &MemEvaluationKeySet{GaloisKeys: map[uint64]*GaloisKey{}}
	eks.RelinearizationKey = relinKey
	for _, k := range galoisKeys {
		eks.GaloisKeys[k.GaloisElement] = k
	}
	return eks
}

func (evk MemEvaluationKeySet) GetGaloisKey(galEl uint64) (gk *GaloisKey, err error) {
	var ok bool
	if gk, ok = evk.GaloisKeys[galEl]; !ok {
		return nil, fmt.Errorf("GaloisKey[%d] is nil", galEl)
	}

	// MINIMAL COMPATIBILITY HOOK: 安全攔截 Fast 密鑰，優雅地回傳 error 阻止誤用
	if gk.Layout == KeyLayoutFast {
		return nil, fmt.Errorf("unsupported: cannot use Fast evaluation key with Standard execution")
	}

	return
}

func (evk MemEvaluationKeySet) GetGaloisKeysList() (galEls []uint64) {
	if evk.GaloisKeys == nil {
		return []uint64{}
	}
	galEls = make([]uint64, len(evk.GaloisKeys))
	var i int
	for galEl := range evk.GaloisKeys {
		galEls[i] = galEl
		i++
	}
	return
}

func (evk MemEvaluationKeySet) GetRelinearizationKey() (rk *RelinearizationKey, err error) {
	if evk.RelinearizationKey != nil {

		// MINIMAL COMPATIBILITY HOOK: 安全攔截 Fast 密鑰，優雅地回傳 error 阻止誤用
		if evk.RelinearizationKey.Layout == KeyLayoutFast {
			return nil, fmt.Errorf("unsupported: cannot use Fast evaluation key with Standard execution")
		}

		return evk.RelinearizationKey, nil
	}

	return nil, fmt.Errorf("RelinearizationKey is nil")
}

func (evk MemEvaluationKeySet) BinarySize() (size int) {
	size++
	if evk.RelinearizationKey != nil {
		size += evk.RelinearizationKey.BinarySize()
	}
	size++
	if evk.GaloisKeys != nil {
		size += evk.GaloisKeys.BinarySize()
	}
	return
}

func (evk MemEvaluationKeySet) WriteTo(w io.Writer) (n int64, err error) {
	switch w := w.(type) {
	case buffer.Writer:
		var inc int64
		if evk.RelinearizationKey != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return inc, err
			}
			n += inc
			if inc, err = evk.RelinearizationKey.WriteTo(w); err != nil {
				return n + inc, err
			}
			n += inc
		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return inc, err
			}
			n += inc
		}

		if evk.GaloisKeys != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return inc, err
			}
			n += inc
			if inc, err = evk.GaloisKeys.WriteTo(w); err != nil {
				return n + inc, err
			}
			n += inc
		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return inc, err
			}
			n += inc
		}
		return n, w.Flush()
	default:
		return evk.WriteTo(bufio.NewWriter(w))
	}
}

func (evk *MemEvaluationKeySet) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:
		var inc int64
		var hasKey uint8
		if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
			return inc, err
		}
		n += inc

		if hasKey == 1 {
			if evk.RelinearizationKey == nil {
				evk.RelinearizationKey = new(RelinearizationKey)
			}
			if inc, err = evk.RelinearizationKey.ReadFrom(r); err != nil {
				return n + inc, err
			}
			n += inc
		}

		if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
			return inc, err
		}
		n += inc

		if hasKey == 1 {
			if evk.GaloisKeys == nil {
				evk.GaloisKeys = structs.Map[uint64, GaloisKey]{}
			}
			if inc, err = evk.GaloisKeys.ReadFrom(r); err != nil {
				return n + inc, err
			}
			n += inc
		}
		return n, nil
	default:
		return evk.ReadFrom(bufio.NewReader(r))
	}
}

func (evk MemEvaluationKeySet) MarshalBinary() (p []byte, err error) {
	buf := buffer.NewBufferSize(evk.BinarySize())
	_, err = evk.WriteTo(buf)
	return buf.Bytes(), err
}

func (evk *MemEvaluationKeySet) UnmarshalBinary(p []byte) (err error) {
	_, err = evk.ReadFrom(buffer.NewBuffer(p))
	return
}
