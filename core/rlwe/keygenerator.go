package rlwe

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/ring/ringqp"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a memory buffer for intermediate values.
type KeyGenerator struct {
	*Encryptor
	pool *BufferPool
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as [EvaluationKey].
func NewKeyGenerator(params ParameterProvider) *KeyGenerator {
	return &KeyGenerator{
		Encryptor: NewEncryptor(params, nil),
		pool:      NewPool(params.GetRLWEParameters().RingQP()),
	}
}

func (kgen KeyGenerator) GenSecretKeyNew() (sk *SecretKey) {
	sk = NewSecretKey(kgen.params)
	kgen.GenSecretKey(sk)
	return
}

func (kgen KeyGenerator) GenSecretKey(sk *SecretKey) {
	kgen.genSecretKeyFromSampler(kgen.xsSampler, sk)
}

func (kgen *KeyGenerator) GenSecretKeyWithHammingWeightNew(hw int) (sk *SecretKey) {
	sk = NewSecretKey(kgen.params)
	kgen.GenSecretKeyWithHammingWeight(hw, sk)
	return
}

func (kgen KeyGenerator) GenSecretKeyWithHammingWeight(hw int, sk *SecretKey) {
	Xs, err := ring.NewSampler(kgen.prng, kgen.params.RingQ(), ring.Ternary{H: hw}, false)
	if err != nil {
		panic(err)
	}
	kgen.genSecretKeyFromSampler(Xs, sk)
}

func (kgen KeyGenerator) genSecretKeyFromSampler(sampler ring.Sampler, sk *SecretKey) {

	ringQP := kgen.params.RingQP().AtLevel(sk.LevelQ(), sk.LevelP())

	sampler.AtLevel(sk.LevelQ()).Read(sk.Value.Q)

	if levelP := sk.LevelP(); levelP > -1 {
		ringQP.ExtendBasisSmallNormAndCenter(sk.Value.Q, levelP, sk.Value.Q, sk.Value.P)
	}

	ringQP.NTT(sk.Value, sk.Value)
	ringQP.MForm(sk.Value, sk.Value)
}

func (kgen KeyGenerator) GenPublicKeyNew(sk *SecretKey) (pk *PublicKey) {
	pk = NewPublicKey(kgen.params)
	kgen.GenPublicKey(sk, pk)
	return
}

func (kgen KeyGenerator) GenPublicKey(sk *SecretKey, pk *PublicKey) {
	pk.Layout = KeyLayoutStandard
	if err := kgen.WithKey(sk).EncryptZero(Element[ringqp.Poly]{
		MetaData: &MetaData{CiphertextMetaData: CiphertextMetaData{IsNTT: true, IsMontgomery: true}},
		Value:    []ringqp.Poly(pk.Value),
	}); err != nil {
		panic(err)
	}
}

func (kgen KeyGenerator) GenKeyPairNew() (sk *SecretKey, pk *PublicKey) {
	sk = kgen.GenSecretKeyNew()
	pk = kgen.GenPublicKeyNew(sk)
	return
}

func (kgen KeyGenerator) GenRelinearizationKeyNew(sk *SecretKey, evkParams ...EvaluationKeyParameters) (rlk *RelinearizationKey) {
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)

	degree := 1
	if compressed {
		degree = 0
	}

	rlk = &RelinearizationKey{EvaluationKey: EvaluationKey{GadgetCiphertext: *NewGadgetCiphertext(kgen.params, degree, levelQ, levelP, BaseTwoDecomposition)}}
	kgen.GenRelinearizationKey(sk, rlk)
	return
}

func (kgen KeyGenerator) GenRelinearizationKey(sk *SecretKey, rlk *RelinearizationKey) {
	sk2 := kgen.pool.AtLevel(rlk.LevelQ()).GetBuffPoly()
	defer kgen.pool.RecycleBuffPoly(sk2)

	sk2.CopyLvl(rlk.LevelQ(), sk.Value.Q)
	kgen.params.RingQ().AtLevel(rlk.LevelQ()).MulCoeffsMontgomery(*sk2, sk.Value.Q, *sk2)
	kgen.genEvaluationKey(*sk2, sk.Value, &rlk.EvaluationKey)
}

func (kgen KeyGenerator) GenGaloisKeyNew(galEl uint64, sk *SecretKey, evkParams ...EvaluationKeyParameters) (gk *GaloisKey) {
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)

	degree := 1
	if compressed {
		degree = 0
	}

	gk = &GaloisKey{
		EvaluationKey: EvaluationKey{GadgetCiphertext: *NewGadgetCiphertext(kgen.params, degree, levelQ, levelP, BaseTwoDecomposition)},
		NthRoot:       kgen.params.GetRLWEParameters().RingQ().NthRoot(),
	}
	kgen.GenGaloisKey(galEl, sk, gk)
	return
}

func (kgen KeyGenerator) GenGaloisKey(galEl uint64, sk *SecretKey, gk *GaloisKey) {

	skIn := sk.Value

	ringQP := kgen.params.RingQP().AtLevel(gk.LevelQ(), gk.LevelP())
	polyQP := kgen.pool.AtLevel(gk.LevelQ(), gk.LevelP())

	skOut := polyQP.GetBuffPolyQP()
	defer polyQP.RecycleBuffPolyQP(skOut)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	galElInv := kgen.params.ModInvGaloisElement(galEl)

	index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), galElInv)

	if err != nil {
		panic(err)
	}

	ringQ.AutomorphismNTTWithIndex(skIn.Q, index, skOut.Q)

	if ringP != nil {
		ringP.AutomorphismNTTWithIndex(skIn.P, index, skOut.P)
	}

	kgen.genEvaluationKey(skIn.Q, *skOut, &gk.EvaluationKey)

	gk.GaloisElement = galEl
	gk.NthRoot = ringQ.NthRoot()
}

func (kgen KeyGenerator) GenGaloisKeys(galEls []uint64, sk *SecretKey, gks []*GaloisKey) {
	if len(galEls) != len(gks) {
		panic(fmt.Errorf("galEls and gks must have the same length"))
	}

	for i, galEl := range galEls {
		if gks[i] == nil {
			gks[i] = kgen.GenGaloisKeyNew(galEl, sk)
		} else {
			kgen.GenGaloisKey(galEl, sk, gks[i])
		}
	}
}

func (kgen KeyGenerator) GenGaloisKeysNew(galEls []uint64, sk *SecretKey, evkParams ...EvaluationKeyParameters) (gks []*GaloisKey) {
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)
	gks = make([]*GaloisKey, len(galEls))
	for i, galEl := range galEls {
		gks[i] = newGaloisKey(kgen.params, levelQ, levelP, BaseTwoDecomposition, compressed)
		kgen.GenGaloisKey(galEl, sk, gks[i])
	}
	return
}

func (kgen KeyGenerator) GenEvaluationKeysForRingSwapNew(skStd, skConjugateInvariant *SecretKey, evkParams ...EvaluationKeyParameters) (stdToci, ciToStd *EvaluationKey) {

	levelQ := utils.Min(skStd.Value.Q.Level(), skConjugateInvariant.Value.Q.Level())

	skCIMappedToStandard := &SecretKey{Value: kgen.params.RingQP().AtLevel(levelQ, kgen.params.MaxLevelP()).NewPoly()}
	kgen.params.RingQ().AtLevel(levelQ).UnfoldConjugateInvariantToStandard(skConjugateInvariant.Value.Q, skCIMappedToStandard.Value.Q)

	if kgen.params.PCount() != 0 {
		buffQ := kgen.pool.GetBuffPoly()
		defer kgen.pool.RecycleBuffPoly(buffQ)
		ExtendBasisSmallNormAndCenterNTTMontgomery(kgen.params.RingQ(), kgen.params.RingP(), skCIMappedToStandard.Value.Q, *buffQ, skCIMappedToStandard.Value.P)
	}

	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)

	stdToci = newEvaluationKey(kgen.params, levelQ, levelP, BaseTwoDecomposition, compressed)
	kgen.GenEvaluationKey(skStd, skCIMappedToStandard, stdToci)

	ciToStd = newEvaluationKey(kgen.params, levelQ, levelP, BaseTwoDecomposition, compressed)
	kgen.GenEvaluationKey(skCIMappedToStandard, skStd, ciToStd)

	return
}

func (kgen KeyGenerator) GenEvaluationKeyNew(skInput, skOutput *SecretKey, evkParams ...EvaluationKeyParameters) (evk *EvaluationKey) {
	levelQ, levelP, BaseTwoDecomposition, compressed := ResolveEvaluationKeyParameters(kgen.params, evkParams)
	evk = newEvaluationKey(kgen.params, levelQ, levelP, BaseTwoDecomposition, compressed)
	kgen.GenEvaluationKey(skInput, skOutput, evk)
	return
}

func (kgen KeyGenerator) GenEvaluationKey(skInput, skOutput *SecretKey, evk *EvaluationKey) {

	ringQ := kgen.params.RingQ()
	ringP := kgen.params.RingP()

	buffSkOut := kgen.pool.GetBuffPolyQP()
	defer kgen.pool.RecycleBuffPolyQP(buffSkOut)
	buffSkIn := kgen.pool.GetBuffPoly()
	defer kgen.pool.RecycleBuffPoly(buffSkIn)

	ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, buffSkOut.Q)

	buffQ := kgen.pool.AtLevel(0).GetBuffPoly()
	defer kgen.pool.RecycleBuffPoly(buffQ)
	if levelP := evk.LevelP(); levelP != -1 {
		ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringP.AtLevel(levelP), buffSkOut.Q, *buffQ, buffSkOut.P)
	}

	ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, *buffSkIn)
	ExtendBasisSmallNormAndCenterNTTMontgomery(ringQ, ringQ.AtLevel(skOutput.Value.Q.Level()), *buffSkIn, *buffQ, *buffSkIn)

	kgen.genEvaluationKey(*buffSkIn, *buffSkOut, evk)
}

func (kgen KeyGenerator) genEvaluationKey(skIn ring.Poly, skOut ringqp.Poly, evk *EvaluationKey) {
	evk.Layout = KeyLayoutStandard

	enc := kgen.WithKey(&SecretKey{Value: skOut})

	if evk.IsCompressed() {
		var seed [32]byte
		if n, err := kgen.prng.Read(seed[:]); n != 32 || err != nil {
			panic(fmt.Errorf("unable to sample evaluation key seed"))
		}
		evk.Seed = &seed

		sampler, err := sampling.NewKeyedPRNG(seed[:])
		if err != nil {
			panic(fmt.Errorf("sampling.NewKeyedPRNG: %w", err))
		}

		enc = enc.withKeyedUniformSampling(sampler)
	}

	for i := 0; i < len(evk.Value); i++ {
		for j := 0; j < len(evk.Value[i]); j++ {
			if err := enc.EncryptZero(Element[ringqp.Poly]{
				MetaData: &MetaData{CiphertextMetaData: CiphertextMetaData{IsNTT: true, IsMontgomery: true}},
				Value:    []ringqp.Poly(evk.Value[i][j]),
			}); err != nil {
				panic(err)
			}
		}
	}

	buffQ := kgen.pool.AtLevel(evk.GadgetCiphertext.LevelQ()).GetBuffPoly()
	defer kgen.pool.RecycleBuffPoly(buffQ)

	if err := AddPolyTimesGadgetVectorToGadgetCiphertext(skIn, []GadgetCiphertext{evk.GadgetCiphertext}, *kgen.params.RingQP(), *buffQ); err != nil {
		panic(err)
	}
}
