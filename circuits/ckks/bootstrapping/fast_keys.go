package bootstrapping

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks/fast"
)

// GenFastEvaluationKeys mirrors Parameters.GenEvaluationKeys while replacing
// only the evaluation-key material with Fast (error, 0) material.
func (p Parameters) GenFastEvaluationKeys(skN1 *rlwe.SecretKey) (btpkeys *EvaluationKeys, skN2 *rlwe.SecretKey, err error) {
	paramsN2 := p.BootstrappingParameters
	kgen := fast.NewKeyGenerator(paramsN2)

	var n1n2, n2n1, realToCmplx, cmplxToReal *rlwe.EvaluationKey
	if p.ResidualParameters.N() != paramsN2.N() {
		skN2 = rlwe.NewKeyGenerator(paramsN2).GenSecretKeyNew()
		if p.ResidualParameters.RingType() == ring.ConjugateInvariant {
			// Match GenEvaluationKeysForRingSwapNew: map the source Q level,
			// but construct the mapped key at the target's full P level.
			levelQ := skN1.LevelQ()
			if skN2.LevelQ() < levelQ {
				levelQ = skN2.LevelQ()
			}
			mapped := paramsN2.RingQP().AtLevel(levelQ, paramsN2.MaxLevelP()).NewPoly()
			paramsN2.RingQ().AtLevel(levelQ).UnfoldConjugateInvariantToStandard(skN1.Value.Q, mapped.Q)
			if paramsN2.PCount() != 0 {
				buff := paramsN2.RingQ().NewPoly()
				rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(paramsN2.RingQ(), paramsN2.RingP(), mapped.Q, buff, mapped.P)
			}
			ciKey := &rlwe.SecretKey{Value: mapped}
			if realToCmplx, err = kgen.GenEvaluationKeyNew(skN2, ciKey); err != nil {
				return
			}
			if cmplxToReal, err = kgen.GenEvaluationKeyNew(ciKey, skN2); err != nil {
				return
			}
		} else {
			if n1n2, err = kgen.GenEvaluationKeyNew(skN1, skN2); err != nil {
				return
			}
			if n2n1, err = kgen.GenEvaluationKeyNew(skN2, skN1); err != nil {
				return
			}
		}
	} else {
		skN2 = rlwe.NewSecretKey(paramsN2)
		buff := paramsN2.RingQ().NewPoly()
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(paramsN2.RingQ(), paramsN2.RingQ(), skN1.Value.Q, buff, skN2.Value.Q)
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(paramsN2.RingQ(), paramsN2.RingP(), skN1.Value.Q, buff, skN2.Value.P)
	}

	var denseToSparse, sparseToDense *rlwe.EvaluationKey
	if p.EphemeralSecretWeight != 0 {
		paramsSparse, _ := rlwe.NewParametersFromLiteral(rlwe.ParametersLiteral{LogN: paramsN2.LogN(), Q: paramsN2.Q()[:1], P: paramsN2.P()[:1]})
		skSparse := rlwe.NewKeyGenerator(paramsSparse).GenSecretKeyWithHammingWeightNew(p.EphemeralSecretWeight)
		kgSparse := fast.NewKeyGenerator(paramsSparse)
		if denseToSparse, err = kgSparse.GenEvaluationKeyNew(skN2, skSparse); err != nil {
			return
		}
		if sparseToDense, err = kgen.GenEvaluationKeyNew(skSparse, skN2); err != nil {
			return
		}
	}

	rlk, err := kgen.GenRelinearizationKeyNew(skN2)
	if err != nil {
		return nil, nil, err
	}
	gks, err := kgen.GenGaloisKeysNew(append(p.GaloisElements(paramsN2), paramsN2.GaloisElementForComplexConjugation()), skN2)
	if err != nil {
		return nil, nil, err
	}
	return &EvaluationKeys{
		EvkN1ToN2: n1n2, EvkN2ToN1: n2n1,
		EvkRealToCmplx: realToCmplx, EvkCmplxToReal: cmplxToReal,
		EvkDenseToSparse: denseToSparse, EvkSparseToDense: sparseToDense,
		MemEvaluationKeySet: rlwe.NewMemEvaluationKeySet(rlk, gks...),
	}, skN2, nil
}

// GenFastBootstrapKeys generates all standard bootstrap keys utilizing Fast architecture (e, 0).
func GenFastBootstrapKeys(
	ckksParams ckks.Parameters,
	sk *rlwe.SecretKey,
	skN1, skN2 *rlwe.SecretKey,
	skReal, skComplex *rlwe.SecretKey,
	skDense, skSparse *rlwe.SecretKey,
	galEls []uint64,
) (
	rlk *rlwe.RelinearizationKey,
	gks []*rlwe.GaloisKey,
	evkN1ToN2 *rlwe.EvaluationKey,
	evkN2ToN1 *rlwe.EvaluationKey,
	evkRealToComplex *rlwe.EvaluationKey,
	evkComplexToReal *rlwe.EvaluationKey,
	evkDenseToSparse *rlwe.EvaluationKey,
	evkSparseToDense *rlwe.EvaluationKey,
	err error,
) {
	kgen := fast.NewKeyGenerator(ckksParams)

	// Relinearization & Galois
	if sk != nil {
		rlk = rlwe.NewRelinearizationKey(ckksParams)
		if err = kgen.GenRelinearizationKey(sk, rlk); err != nil {
			return
		}
		if len(galEls) > 0 {
			if gks, err = kgen.GenGaloisKeysNew(galEls, sk); err != nil {
				return
			}
		}
	}

	// Ring Degree Switches (N1 <-> N2)
	if skN1 != nil && skN2 != nil {
		if evkN1ToN2, err = kgen.GenEvaluationKeyNew(skN1, skN2); err != nil {
			return
		}
		if evkN2ToN1, err = kgen.GenEvaluationKeyNew(skN2, skN1); err != nil {
			return
		}
	}

	// Ring Type Switches (Real <-> Complex)
	if skReal != nil && skComplex != nil {
		if evkRealToComplex, err = kgen.GenEvaluationKeyNew(skReal, skComplex); err != nil {
			return
		}
		if evkComplexToReal, err = kgen.GenEvaluationKeyNew(skComplex, skReal); err != nil {
			return
		}
	}

	// Slot Structure Switches (Dense <-> Sparse)
	if skDense != nil && skSparse != nil {
		if evkDenseToSparse, err = kgen.GenEvaluationKeyNew(skDense, skSparse); err != nil {
			return
		}
		if evkSparseToDense, err = kgen.GenEvaluationKeyNew(skSparse, skDense); err != nil {
			return
		}
	}

	return
}
