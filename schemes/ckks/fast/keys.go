package fast

import (
	"fmt"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// KeyGenerator wraps the standard RLWE key generator to produce Fast keys.
type KeyGenerator struct {
	*rlwe.KeyGenerator
	params rlwe.ParameterProvider
}

// NewKeyGenerator instantiates a KeyGenerator that exclusively produces Fast keys.
func NewKeyGenerator(params rlwe.ParameterProvider) *KeyGenerator {
	return &KeyGenerator{
		KeyGenerator: rlwe.NewKeyGenerator(params),
		params:       params,
	}
}

// getZeroKey returns a SecretKey with all elements forced to exactly zero
func getZeroKey(sk *rlwe.SecretKey) *rlwe.SecretKey {
	skZ := sk.CopyNew()
	skZ.Value.Q.Zero()
	if skZ.Value.P.Coeffs != nil {
		skZ.Value.P.Zero()
	}
	return skZ
}

func (kgen *KeyGenerator) GenPublicKeyNew(sk *rlwe.SecretKey) *rlwe.PublicKey {
	pk := rlwe.NewPublicKey(kgen.params)
	kgen.GenPublicKey(sk, pk)
	return pk
}

// GenPublicKey creates a Fast PublicKey (e, 0)
func (kgen *KeyGenerator) GenPublicKey(sk *rlwe.SecretKey, pk *rlwe.PublicKey) {
	skZero := getZeroKey(sk)
	kgen.KeyGenerator.GenPublicKey(skZero, pk)
	pk.Layout = rlwe.KeyLayoutFast
	if len(pk.Value) > 1 {
		pk.Value[1].Q.Zero()
		if pk.Value[1].P.Coeffs != nil {
			pk.Value[1].P.Zero()
		}
	}
}

func (kgen *KeyGenerator) fastMaterial(evk *rlwe.EvaluationKey, skIn, skOut *rlwe.SecretKey) error {
	if evk.IsCompressed() {
		return fmt.Errorf("unsupported: Fast EvaluationKeys cannot be compressed")
	}
	zeroIn, zeroOut := getZeroKey(skIn), getZeroKey(skOut)
	kgen.KeyGenerator.GenEvaluationKey(zeroIn, zeroOut, evk)
	evk.Layout = rlwe.KeyLayoutFast
	evk.GadgetCiphertext.Layout = rlwe.KeyLayoutFast
	for i := range evk.Value {
		for j := range evk.Value[i] {
			evk.Value[i][j][1].Q.Zero()
			if evk.Value[i][j][1].P.Coeffs != nil {
				evk.Value[i][j][1].P.Zero()
			}
		}
	}
	return nil
}

func (kgen *KeyGenerator) GenEvaluationKeyNew(skInput, skOutput *rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (*rlwe.EvaluationKey, error) {
	evk := rlwe.NewEvaluationKey(kgen.params, evkParams...)
	return evk, kgen.GenEvaluationKey(skInput, skOutput, evk)
}

// GenEvaluationKey intercepts EvaluationKey generation, forcing inputs to 0 to extract native error.
func (kgen *KeyGenerator) GenEvaluationKey(skInput, skOutput *rlwe.SecretKey, evk *rlwe.EvaluationKey) error {
	return kgen.fastMaterial(evk, skInput, skOutput)
}

func (kgen *KeyGenerator) GenRelinearizationKeyNew(sk *rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (*rlwe.RelinearizationKey, error) {
	rlk := rlwe.NewRelinearizationKey(kgen.params, evkParams...)
	return rlk, kgen.GenRelinearizationKey(sk, rlk)
}

func (kgen *KeyGenerator) GenRelinearizationKey(sk *rlwe.SecretKey, rlk *rlwe.RelinearizationKey) error {
	// Preserve the standard sk^2 construction before replacing only key material.
	std := rlwe.NewKeyGenerator(kgen.params)
	std.GenRelinearizationKey(sk, rlk)
	return kgen.fastMaterial(&rlk.EvaluationKey, sk, sk)
}

func (kgen *KeyGenerator) GenGaloisKeyNew(galEl uint64, sk *rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (*rlwe.GaloisKey, error) {
	gk := rlwe.NewGaloisKey(kgen.params, evkParams...)
	return gk, kgen.GenGaloisKey(galEl, sk, gk)
}

func (kgen *KeyGenerator) GenGaloisKey(galEl uint64, sk *rlwe.SecretKey, gk *rlwe.GaloisKey) error {
	// Run the standard path first to preserve inverse-Galois and metadata semantics.
	std := rlwe.NewKeyGenerator(kgen.params)
	std.GenGaloisKey(galEl, sk, gk)
	return kgen.GenEvaluationKey(sk, sk, &gk.EvaluationKey)
}

func (kgen *KeyGenerator) GenGaloisKeysNew(galEls []uint64, sk *rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) ([]*rlwe.GaloisKey, error) {
	gks := make([]*rlwe.GaloisKey, len(galEls))
	for i, galEl := range galEls {
		var err error
		if gks[i], err = kgen.GenGaloisKeyNew(galEl, sk, evkParams...); err != nil {
			return nil, err
		}
	}
	return gks, nil
}
