# Fast Bootstrap Coverage and Performance Gap Audit

Status: read-only source audit of the `fast-ckks` branch at `410d8743`

Authority: subordinate to `docs/FAST_CKKS_SPEC.md`

Stage: Stage A — complete Fast CKKS/Bootstrap execution with maximum practical speed

## 1. Executive conclusion

The branch has useful Fast primitives, but the public Bootstrap execution path does not use them. `bootstrapping.NewEvaluator` constructs a Standard `*ckks.Evaluator`, then builds the DFT and Mod1 evaluators around that Standard evaluator. Consequently, the live path performs full-Q arithmetic, QP decomposition, GadgetProduct, key switching, ModDown, and Standard Rescale throughout.

The shortest path is not to rewrite Bootstrap. It is to supply the missing Fast evaluator contracts in dependency order and wire the existing orchestration to them:

1. define and implement Fast q0/q1 level-transition/Rescale semantics;
2. provide an NTT-compatible Fast evaluator surface for polynomial circuits;
3. adapt DFT to the existing Fast LinearTransform/automorphism primitives;
4. define the deliberate q0/q1 ModUp boundary;
5. wire packing/ring-degree conversion and the Bootstrap orchestrator.

The top architectural conflict is that Standard `ScaleDown` reduces to level 0 and Standard `ModUp` reconstructs all Q limbs from q0, while established Fast primitives require q0 and q1 and treat q2...qL as dormant. This must be resolved explicitly; full redistribution is not an acceptable hidden workaround.

## 2. Audit scope and notation

This audit traced only the current Bootstrap call graph and the directly reached CKKS, DFT, polynomial, ring, and RLWE operations. Completed Fast primitives were inspected only to identify wiring and contract gaps.

- A ciphertext at level `L` stores `L+1` Q limbs, `q0...qL`.
- Fast authority is q0/q1; q2...qL are dormant.
- `c0`, `c1`, and `c2` are ciphertext polynomial components, not RNS limbs.
- Cost classifications are structural; no benchmark percentages are claimed.

## 3. Actual public execution graph

### 3.1 Standard-ring public path

```text
circuits/ckks/bootstrapping/evaluator.go
Evaluator.Bootstrap
└── Evaluator.BootstrapMany
    ├── Evaluator.PackAndSwitchN1ToN2
    │   ├── Evaluator.pack (N1, when N1 != N2)
    │   │   └── full-Q MulCoeffsMontgomeryThenAdd
    │   ├── Evaluator.switchRingDegreeN1ToN2New (when N1 != N2)
    │   │   └── ckks.Evaluator.ApplyEvaluationKey
    │   │       └── rlwe.Evaluator.ApplyEvaluationKey
    │   │           └── GadgetProduct -> QP work -> ModDown
    │   └── Evaluator.pack (N2)
    │       └── full-Q MulCoeffsMontgomeryThenAdd
    ├── Evaluator.Evaluate (once per packed ciphertext)
    │   ├── Evaluator.bootstrap
    │   │   ├── Evaluator.ScaleDown
    │   │   │   ├── ckks.Evaluator.Mul (scalar)
    │   │   │   └── ckks.Evaluator.RescaleTo (when level != 0)
    │   │   ├── Evaluator.ModUp
    │   │   │   ├── ApplyEvaluationKey dense->sparse (optional)
    │   │   │   ├── full-Q INTT
    │   │   │   ├── q0-centered materialization into q1...qL
    │   │   │   ├── QP decomposition/NTT/GadgetProductHoisted (optional sparse->dense)
    │   │   │   ├── full-Q NTT and scalar multiplication
    │   │   │   └── rlwe.Evaluator.Trace
    │   │   │       └── repeated Automorphism -> GadgetProduct/ModDown + full-Q Add
    │   │   ├── Evaluator.CoeffsToSlots
    │   │   │   └── dft.Evaluator.CoeffsToSlotsNew
    │   │   │       ├── dft.Evaluator.dft
    │   │   │       │   ├── common/lintrans.Evaluator.Evaluate
    │   │   │       │   │   ├── DecomposeNTT into QP
    │   │   │       │   │   ├── hoisted/BSGS Automorphism
    │   │   │       │   │   ├── GadgetProductHoistedLazy
    │   │   │       │   │   ├── full-QP plaintext products
    │   │   │       │   │   └── ModDownQPtoQNTT
    │   │   │       │   └── ckks.Evaluator.Rescale per factorization group
    │   │   │       ├── Conjugate, Add/Sub, scalar Mul
    │   │   │       └── optional Rotate for imaginary repacking
    │   │   ├── Evaluator.EvalMod (real and, for dense packing, imaginary)
    │   │   │   └── mod1.Evaluator.EvaluateNew
    │   │   │       ├── polynomial.Evaluator.Evaluate
    │   │   │       │   └── common/polynomial Paterson-Stockmeyer evaluation
    │   │   │       │       ├── PowerBasis.GenPower
    │   │   │       │       ├── Mul/MulRelin
    │   │   │       │       ├── Relinearize
    │   │   │       │       ├── repeated Rescale
    │   │   │       │       ├── scalar/plaintext MulThenAdd
    │   │   │       │       └── Add/Sub
    │   │   │       └── DoubleAngle loop: MulRelin, Add, scalar Add, Rescale
    │   │   └── Evaluator.SlotsToCoeffs
    │   │       └── dft.Evaluator.SlotsToCoeffsNew
    │   │           ├── optional scalar Mul + Add
    │   │           └── dft.Evaluator.dft (same LT/QP/Rescale chain)
    │   └── optional iterative-precision loop
    │       └── repeated bootstrap plus Standard Mul/Sub/Rescale
    └── Evaluator.UnpackAndSwitchN2ToN1
        ├── Evaluator.unpack (N2): full-Q monomial multiplication
        ├── switchRingDegreeN2ToN1New (when N1 != N2)
        │   └── ApplyEvaluationKey -> GadgetProduct/ModDown
        └── Evaluator.unpack (N1, when N1 != N2)
```

### 3.2 Configuration-dependent branches

| Condition | Current branch |
|---|---|
| Residual ring is `ring.Standard` | Uses packing, optional N1/N2 switching, `Evaluate`, then unpacking |
| Residual ring is `ring.ConjugateInvariant` | `EvaluateConjugateInvariant`: Standard RealToComplex, optional imaginary packing, `Evaluate`, then Standard ComplexToReal |
| `N1 != N2` | Standard evaluation-key ring switch before and after Bootstrap |
| Ephemeral secret enabled | ModUp uses dense/sparse evaluation keys, QP decomposition, and hoisted GadgetProduct |
| Dense slot packing | CoeffsToSlots emits real and imaginary ciphertexts; EvalMod runs twice |
| DFT matrix has `N1 != 0` BSGS split | Standard common LinearTransform uses BSGS/double hoisting; otherwise it uses the naive hoisted path |
| Iteration parameters enabled | `Evaluate` invokes additional Bootstrap rounds and Standard correction arithmetic |

`parameters_literal.go` defines `ModUpThenEncode`, `DecodeThenModUp`, and `Custom`. However, the public `Evaluator.bootstrap` method in this branch always executes `ScaleDown -> ModUp -> CoeffsToSlots -> EvalMod -> SlotsToCoeffs`. `CircuitOrder` changes constructor consistency checks, but no dispatch on `CircuitOrder` exists in `bootstrap`. Therefore `DecodeThenModUp` is not a distinct public `Bootstrap` execution graph in the audited source; users would need manual/custom stage orchestration. This discrepancy is `UNKNOWN_REQUIRES_TARGETED_AUDIT` before relying on that order.

## 4. Source-node inventory and classification

Each operation has exactly one Fast coverage classification.

| Source node | Caller and important downstream work | Classification | Why |
|---|---|---|---|
| `bootstrapping.Evaluator.Bootstrap/BootstrapMany` in `evaluator.go` | Public entry; pack, Evaluate, unpack | `MISSING_FAST_IMPLEMENTATION` | No Fast orchestrator or evaluator selection exists |
| `bootstrapping.NewEvaluator` in `evaluator.go` | Constructs `ckks.NewEvaluator`, `dft.NewEvaluator`, and `mod1.NewEvaluator` | `STANDARD_HIGH_COST` | This binds every circuit stage to Standard execution |
| `PackAndSwitchN1ToN2.pack` / `unpack` | Full-Q monomial multiply/add | `MISSING_FAST_IMPLEMENTATION` | No q0/q1 packing wrapper; dormant limbs are processed at active level |
| `switchRingDegreeN1ToN2New` / `switchRingDegreeN2ToN1New` | Standard `ApplyEvaluationKey` | `FAST_EXISTS_BUT_NOT_WIRED` | `FastN1ToN2/FastN2ToN1` exist, but Bootstrap calls Standard key switching |
| `ScaleDown` level dropping | `Resize` only | `STANDARD_BUT_CHEAP` | Metadata/storage truncation is cheap, but dropping to level 0 conflicts with q0/q1 authority |
| `ScaleDown` scalar multiply | Standard CKKS scalar Mul over all active Q limbs | `MISSING_FAST_IMPLEMENTATION` | A q0/q1 scalar path is absent from the Fast evaluator surface |
| `ScaleDown` `RescaleTo` | Full-Q rounded division and level consumption | `MISSING_FAST_IMPLEMENTATION` | Fast Rescale is absent and is a Bootstrap blocker |
| `ModUp` | INTT, q0-to-all-Q materialization, NTT, optional QP key switch, Trace | `BOUNDARY_ONLY` | Modulus raising is an explicit representation boundary with unresolved q0/q1 semantics |
| `ModUp` dense/sparse switches | ApplyEvaluationKey and GadgetProductHoisted | `STANDARD_HIGH_COST` | Default-like ephemeral-secret configurations activate QP security work |
| `Trace` | Scalar multiply, repeated Standard automorphisms and adds | `FAST_EXISTS_BUT_NOT_WIRED` | Fast automorphism and q0/q1 Add exist; a Fast Trace adapter does not |
| `dft.Evaluator.dft` | Repeated common LinearTransform then Rescale | `STANDARD_HIGH_COST` | Entire DFT remains full-Q/QP and repeats by factorization depth |
| common `LinearTransform.EvaluateMany` | DecomposeNTT, hoisted/BSGS rotations, QP products, ModDown | `FAST_EXISTS_BUT_NOT_WIRED` | Fast LinearTransform exists but its level/Montgomery/API contracts do not match this caller |
| DFT Conjugate/Rotate | Standard automorphism/key switching | `FAST_EXISTS_BUT_NOT_WIRED` | Fast Standard-ring automorphism exists, but DFT calls Standard CKKS methods |
| DFT Add/Sub | Standard full-Q arithmetic | `FAST_EXISTS_BUT_NOT_WIRED` | FastAdd/FastSub exist but are not exposed through DFT's evaluator |
| DFT scalar Mul | Standard full-Q scalar operation | `MISSING_FAST_IMPLEMENTATION` | Current Fast Mul only accepts two ciphertexts |
| DFT Rescale | Standard full-Q rounded division | `MISSING_FAST_IMPLEMENTATION` | Required after each DFT factorization group |
| `mod1.Evaluator.EvaluateNew` | Polynomial evaluation and DoubleAngle | `STANDARD_HIGH_COST` | It is constructed with Standard CKKS and polynomial evaluators |
| Polynomial ciphertext Mul/MulRelin | Power basis and Paterson-Stockmeyer steps | `FAST_EXISTS_BUT_NOT_WIRED` | Fast ciphertext Mul exists, but only for coefficient/non-Montgomery inputs and not the generic evaluator interface |
| Polynomial Relinearize | Standard GadgetProduct | `FAST_EXISTS_BUT_NOT_WIRED` | Current zero-secret truncation exists but is not wired as `Relinearize` |
| Polynomial Rescale | Repeated full-Q rounded division | `MISSING_FAST_IMPLEMENTATION` | No Fast equivalent; repeatedly blocks EvalMod |
| Polynomial scalar/plaintext Mul/MulThenAdd | Coefficients and correction terms | `MISSING_FAST_IMPLEMENTATION` | General Fast evaluator forms are absent |
| Polynomial Add/Sub | Full-Q Standard arithmetic | `FAST_EXISTS_BUT_NOT_WIRED` | FastAdd/FastSub exist but not behind `schemes.Evaluator` |
| `GenFastEvaluationKeys` / `GenFastBootstrapKeys` in `fast_keys.go` | Generates Fast-layout Bootstrap key material | `ALREADY_FAST` | Key support exists, but Standard consumers explicitly reject Fast layouts |
| Final scale assignment and `DropLevel` | Metadata-only operations | `STANDARD_BUT_CHEAP` | No heavy arithmetic |
| Conjugate-invariant ring switching | Standard DomainSwitcher and evaluation keys | `UNKNOWN_REQUIRES_TARGETED_AUDIT` | Current Fast automorphism/ring-degree code supports Standard rings only |

## 5. Hidden Standard fallbacks

### HIDDEN_STANDARD_FALLBACK — Bootstrap evaluator construction

- Path: `circuits/ckks/bootstrapping/evaluator.go`, `NewEvaluator`.
- Evidence: `eval.Evaluator = ckks.NewEvaluator(params, evk)`; both `dft.NewEvaluator` and `mod1.NewEvaluator` receive this Standard evaluator.
- Effect: every reached DFT, Mod1, polynomial, rotation, relinearization, Rescale, and key-switch operation uses Standard full-Q/QP behavior.

### HIDDEN_STANDARD_FALLBACK — DFT LinearTransform

- Path: `circuits/ckks/dft/dft.go:dft` -> `circuits/common/lintrans/lintrans_evaluator.go:EvaluateMany`.
- Effect: q0/q1 Fast LinearTransform is bypassed; Standard QP decomposition, GadgetProductHoistedLazy, BSGS/hoisting, and ModDown execute for each factorized matrix.

### HIDDEN_STANDARD_FALLBACK — EvalMod polynomial evaluator

- Path: `bootstrapping.NewEvaluator` -> `mod1.NewEvaluator` -> `polynomial.NewEvaluator` -> `circuits/common/polynomial`.
- Effect: current Fast Mul and truncation are bypassed; Standard MulRelin/Relinearize/GadgetProduct and full-Q Rescale execute repeatedly.

### HIDDEN_STANDARD_FALLBACK — ring switching and ModUp key switches

- Paths: `switchRingDegreeN1ToN2New`, `switchRingDegreeN2ToN1New`, and optional dense/sparse branches in `ModUp`.
- Effect: Standard `ApplyEvaluationKey` and GadgetProduct are selected instead of existing q0/q1 ring conversion or a Fast key-switch contract.

This fallback is also a functional incompatibility with Fast keys: `core/rlwe/evaluator_evaluationkey.go:ApplyEvaluationKey` and `core/rlwe/evaluator_gadget_product.go:GadgetProduct/GadgetProductLazy` explicitly reject `KeyLayoutFast`. Fast key generation alone cannot make Standard Bootstrap run Fast.

## 6. Remaining full-RNS work

| Location | Operation and limbs | Are q2...qL necessary for current Stage A? |
|---|---|---|
| `BootstrapMany` packing/unpacking | `RingQ.AtLevel(level).MulCoeffsMontgomery*` on all `level+1` Q limbs | They appear to be Standard representation maintenance; q0/q1 monomial multiplication is sufficient in the existing Fast model, subject to adapter validation |
| `ScaleDown` scalar Mul | All Q limbs at the input level | Dormant-limb work appears unnecessary; q0/q1 scalar arithmetic is missing |
| `ckks.Evaluator.Rescale/RescaleTo` | For each component, consumes highest qL and updates every remaining q0...q(L-1); one rescale may consume multiple primes | Standard exact rounded semantics need the dropped limb. With q2...qL dormant, the required Fast semantic is unresolved; full-Q work cannot simply be retained |
| `ModUp` | Full INTT at input level, writes every q1...qL from centered q0, then full-Q NTT and scalar Mul | Full materialization is a Standard Bootstrap invariant, not an established Fast requirement. The current q0-only source conflicts with q0/q1 authority |
| `Trace` | Full-Q scalar multiply, automorphisms, and additions at max Bootstrap level | Fast q0/q1 primitives appear sufficient under current zero-secret semantics; wiring is absent |
| DFT common LinearTransform | Every active Q limb plus every configured P limb; repeated for all factor matrices | q2...qL appear unnecessary for current Fast authority; Standard code maintains them because its ciphertext state is full-RNS |
| DFT Rescale | All remaining Q limbs after each factorization group | Same unresolved Fast Rescale issue; invoked `len(Matrix.Levels)` times per DFT |
| EvalMod polynomial arithmetic | Full-Q Mul/Add/MulThenAdd at each current level | Dormant-limb work appears unnecessary except where Standard Rescale depends on dropped limbs |
| EvalMod relinearization | Full-Q output plus QP gadget path | Candidate for direct degree truncation in current zero-secret mode |
| CI domain switching | Full active Q and evaluation-key work | Uncertain because existing Fast code does not support the conjugate-invariant mapping |

## 7. QP/P path audit

| Critical-path use | Source | Classification | Rationale |
|---|---|---|---|
| DFT rotations and LinearTransform | `circuits/common/lintrans/lintrans_evaluator.go` | `STANDARD_SECURITY_PATH_CANDIDATE_FOR_ELISION` | Existing Fast automorphism and LinearTransform avoid keys, QP, decomposition, and ModDown in current zero-secret mode |
| EvalMod MulRelin/Relinearize | `schemes/ckks/evaluator.go:mulRelin`; `core/rlwe/evaluator_evaluationkey.go:Relinearize` | `STANDARD_SECURITY_PATH_CANDIDATE_FOR_ELISION` | Existing Fast multiplication discards c2 under the current mode |
| ModUp dense->sparse and sparse->dense | `bootstrapping/evaluator.go:ModUp` | `STANDARD_SECURITY_PATH_CANDIDATE_FOR_ELISION` | Security-oriented ephemeral-secret switching is not a Stage-A objective; exact removal must be configuration-aware |
| Trace/Conjugate/Rotate | `core/rlwe/evaluator_automorphism.go` | `STANDARD_SECURITY_PATH_CANDIDATE_FOR_ELISION` | Standard automorphism invokes GadgetProduct; current Fast automorphism does not |
| N1/N2 ring-degree switch | `core/rlwe/evaluator_evaluationkey.go:ApplyEvaluationKey` | `STANDARD_SECURITY_PATH_CANDIDATE_FOR_ELISION` | Existing Fast q0/q1 ring-degree conversion is designed to avoid the evaluation key |
| Conjugate-invariant Real/Complex switch | CKKS DomainSwitcher reached by `EvaluateConjugateInvariant` | `UNCERTAIN` | Current Fast Standard-ring-only primitives do not establish equivalent semantics |
| Future nonzero-secret/noise experiments | Fast specification extension points | `UNCERTAIN` | P, gadget terms, and key switching must remain structurally recoverable for Stage B audit |

No QP/P operation is proven `REQUIRED_BY_CURRENT_FAST_STAGE_A` by the current zero-secret architecture. ModUp's sparse-key branch and CI switching remain audit-sensitive boundaries; this document does not authorize their removal.

## 8. Rescale audit

### Standard implementation

- Entry points: `schemes/ckks/evaluator.go:Rescale` and `RescaleTo`.
- Expected input: ciphertext metadata present; the arithmetic path calls `ring.DivRoundByLastModulusManyNTT`, whose contract requires NTT-domain input.
- Scale update: divide `Scale` by each consumed trailing Q prime.
- Level update: decrease by `LevelsConsumedPerRescaling()` for `Rescale`, or by the number selected against `minScale` for `RescaleTo`.
- Limb access: at level `L`, `DivRoundByLastModulusNTT` reads qL, inverse-transforms that limb, and updates each lower q0...q(L-1). Multiple-prime rescale performs full INTT, sequential divisions, and full NTT on the remaining limbs.
- Allocation: one pool polynomial per ciphertext Rescale plus ring-internal scratch; multi-prime paths allocate another ring buffer.

### Bootstrap frequency

- `ScaleDown` invokes `RescaleTo` when its input remains above level 0.
- Each DFT invokes `Rescale` once per outer factorization group (`len(Matrix.Levels)`), for both CoeffsToSlots and SlotsToCoeffs. Repository defaults define four CoeffsToSlots groups and three SlotsToCoeffs groups, but active parameters are configurable.
- Power-basis generation and Paterson-Stockmeyer combination invoke Rescale repeatedly according to polynomial depth.
- Mod1's DoubleAngle loop invokes Rescale once per iteration (default literal: three).
- Optional precision iterations can invoke another Rescale and additional complete Bootstrap rounds.

### Fast conclusion

Metadata and level changes are easy to preserve, but Standard rounded division fundamentally reads the prime being dropped. If q2...qL are dormant, their residues cannot be used as Standard rounding inputs. A Fast q0/q1 rule must specify how authoritative residues and scale evolve without silently reviving full-RNS state.

Classification: `BOOTSTRAP_BLOCKER`. Priority: `HIGH`.

## 9. ModUp / modulus-raising audit

`bootstrapping.Evaluator.ModUp` begins after `ScaleDown`. Standard `ScaleDown` aims for level 0, so the normal ModUp source is q0, not q0/q1.

The current path:

1. optionally applies a dense-to-sparse evaluation key;
2. applies INTT to every current ciphertext component over all active input limbs;
3. resizes storage to Bootstrap `MaxLevel`;
4. centers each q0 coefficient and reduces it into q1...qL for c0 and, without ephemeral switching, c1;
5. with ephemeral switching, constructs QP decomposition buffers, performs full-Q/full-P NTT and scalar multiplication, then calls `GadgetProductHoisted`;
6. otherwise performs full-Q NTT and scalar multiplication directly;
7. calls `Trace`, which performs repeated Standard key-switched automorphisms.

This is full redistribution on the Standard hot path. It conflicts with Fast's q0/q1-authoritative model in two ways: it starts from q0-only state, and it materializes every dormant Q limb before all later work.

The existing `ReconstructQ0Q1`/`Redistribute` helper is a reference/debug oracle. Its q0/q1 reconstruction may help define an explicit boundary, but its big.Int reconstruction and full redistribution are not suitable for production Fast execution.

Unresolved requirements:

- whether Fast ScaleDown should stop at level 1 to retain q0/q1;
- whether Stage-A ModUp should be metadata-only, q0/q1 normalization, or another bounded transformation;
- how the q0/q1 value and scale must be adjusted before CoeffsToSlots/EvalMod;
- which sparse-packing Trace semantics are required without Standard key switching.

Classification: `BOUNDARY_ONLY`. Performance priority: `HIGH`. Mathematical/representation decision required before implementation.

## 10. CoeffsToSlots and SlotsToCoeffs

Both transforms call `dft.Evaluator.dft`. For every matrix in each factorization group, `dft` calls the common LinearTransform evaluator; after the group, it calls Standard Rescale. The common evaluator chooses naive hoisting or BSGS based on matrix `N1`, decomposes c1 into QP once, performs key-switched automorphisms, multiplies QP plaintext diagonals, accumulates, and ModDowns QP to Q.

Additional CoeffsToSlots operations include Conjugate, Add/Sub, scalar multiplication by `-i`, and an optional Rotate for imaginary repacking. SlotsToCoeffs may multiply the imaginary ciphertext by `i` and add it before its DFT.

### Can current Fast LinearTransform replace a Standard call directly?

No, not from the current DFT caller.

Exact contract mismatches:

- **Evaluator API:** DFT owns a common LinearTransform evaluator built around `schemes.Evaluator`; `fast.Evaluator` does not implement that interface.
- **Level:** Fast LinearTransform requires `matrix.LevelQ == ctIn.Level()` and equal input/output levels. DFT consumes levels between factorization groups, while matrix construction/caller selection permits evaluation at the current lower level.
- **Montgomery representation:** Fast LinearTransform requires both ciphertext and matrix to be Montgomery. Normal CKKS ciphertexts are NTT by default but not marked Montgomery; encoded matrices are NTT/Montgomery.
- **Rescale:** Fast LinearTransform performs no level transition; DFT requires one Rescale per factorization group.
- **BSGS/hoisting:** Fast LinearTransform loops over all diagonals with direct Fast automorphisms. It has no BSGS/hoisted interface. This is functionally plausible for Stage A but may be slower for large diagonal sets until measured.
- **Ring type:** Fast LinearTransform supports `ring.Standard` only.
- **Scratch:** it allocates six q0/q1 polynomials per call, and each Fast automorphism allocates temporary N-length slices.

The mismatch is not missing Fast diagonal arithmetic; it is a missing DFT adapter plus level/domain/representation contracts.

## 11. EvalMod / Mod1 audit

The active Mod1 path evaluates a Chebyshev polynomial using the generic Paterson-Stockmeyer evaluator. It generates a power basis with ciphertext Mul/MulRelin, conditionally keeps degree-2 intermediates, relinearizes them, rescales powers and giant steps, multiplies powers by scalar/plaintext coefficients, and accumulates with Add/Sub. The DoubleAngle loop then repeats `MulRelin(res,res)`, additions, and `Rescale`.

The existing Fast Mul does not satisfy these calls directly:

- it accepts only two degree-1 ciphertext pointers, not the `rlwe.Operand` forms required by `schemes.Evaluator`;
- it rejects NTT or Montgomery inputs, while Bootstrap polynomial evaluation remains NTT-domain;
- it always returns a degree-1 result through current zero-secret truncation, while generic lazy power-basis code sometimes expects temporary degree 2 and calls `Relinearize` later;
- Fast `Relinearize` and `Rescale` methods with the generic signatures do not exist;
- Fast scalar/plaintext `Mul`, `MulThenAdd`, and `Add/Sub` adapter methods are absent.

Current Fast multiplication arithmetic is reusable, but its representation and interface are blockers. EvalMod is `STANDARD_HIGH_COST`; Fast evaluator coverage for it is a `BLOCKER` with `HIGH` priority.

## 12. Ring-degree switching

When `ResidualParameters.N() != BootstrappingParameters.N()`, packing calls Standard `ApplyEvaluationKey` in both directions. This performs ring mapping plus evaluation-key GadgetProduct and, for a down-switch, a temporary ciphertext.

`FastN1ToN2` and `FastN2ToN1` already provide q0/q1-only ring mapping without keys, CRT, or redistribution. They are **partially sufficient and need an adapter**:

- Bootstrap must call them instead of `ApplyEvaluationKey` under the current Fast mode;
- input/output rings must expose matching q0/q1 moduli and equal active levels;
- NTT inputs trigger partial INTT -> coefficient mapping -> partial NTT, creating `DOMAIN_CONVERSION_DEBT`;
- each NTT conversion allocates two polynomials per ciphertext component;
- only Standard rings are supported.

This is not a missing mathematical primitive for the Standard-ring N1/N2 case.

## 13. Domain-conversion map

| Boundary/stage | Current expected representation | Fast issue |
|---|---|---|
| Bootstrap input and pack | Normally NTT; ciphertext is generally non-Montgomery | Pack uses full-Q NTT-domain monomial products |
| ScaleDown | NTT through scalar Mul and RescaleTo | Current Fast Mul cannot accept NTT inputs |
| ModUp entry | NTT | Standard path performs full INTT |
| ModUp middle | Coefficient, non-Montgomery | q0 is centered and redistributed to all Q/P limbs |
| ModUp exit and Trace | NTT, non-Montgomery ciphertext | Standard performs full-Q NTT; Trace rotations key-switch in QP |
| DFT matrices | NTT and Montgomery plaintext diagonals | Compatible with Fast diagonal storage, but Fast LT also requires the ciphertext to be marked Montgomery |
| DFT ciphertext | NTT, normally non-Montgomery | Direct Fast LT rejects it |
| EvalMod | NTT across polynomial arithmetic and Rescale | Current Fast ciphertext Mul rejects it and internally performs coefficient->NTT->coefficient work |
| Ring-degree switch | Preserves input domain | Existing Fast NTT switch internally performs partial INTT and NTT |

`DOMAIN_CONVERSION_DEBT`:

- Adapting current Fast Mul to Bootstrap by wrapping every call in NTT -> INTT -> Fast Mul -> NTT would multiply conversion cost across the polynomial evaluator and is not acceptable as the production design.
- `FastMulQ01Authoritative` currently invokes partial NTT and partial INTT for each polynomial product. `fast.Evaluator.Mul` invokes it three times per ciphertext multiplication.
- NTT ring-degree conversion performs partial INTT/NTT for both c0 and c1.
- ModUp's full INTT/full NTT pair is a major Standard boundary cost.

The preferred Stage-A direction is an NTT-resident q0/q1 evaluator across DFT and EvalMod, with deliberate conversions only at true boundaries.

## 14. Allocation map

`ALLOCATION_DEBT` identified from source:

- `fast.Evaluator.Mul` allocates three q0/q1 result temporaries; its three polynomial products each allocate four more q0/q1 polynomials in `fastMulQ01Core` (15 polynomial allocations per ciphertext multiplication before output allocation).
- `fast.Evaluator.Automorphism` allocates two q0/q1 polynomials; `FastAutomorphism` additionally allocates one N-length `[]uint64` per limb per polynomial call.
- `fast.Evaluator.LinearTransform` allocates six q0/q1 polynomials per call and invokes two allocation-heavy automorphisms for every diagonal.
- Fast NTT ring-degree conversion allocates input/output coefficient polynomials for each of c0 and c1.
- Standard common LinearTransform allocates QP decomposition, several QP scratch polynomials, output buffers, and pre-rotated extended ciphertexts; pools reduce churn but the structures remain large.
- DFT allocates output ciphertexts and, for split real/imag handling, copies or borrows temporary ciphertexts.
- Polynomial `PowerBasis` stores multiple ciphertext powers; Paterson-Stockmeyer allocates baby-step slices and result ciphertexts.
- Bootstrap packing/unpacking copies ciphertexts and grows output slices; optional iterative precision invokes additional full Bootstrap allocations.

These are secondary to coverage and structural elision. Do not optimize them before the Fast path is complete enough to benchmark end-to-end.

## 15. Consolidated implementation-order gap table

| Order | Bootstrap stage | Operation | Current path | Fast status | Main expensive work | Blocker? | Recommended next task |
| ----: | --------------- | --------- | ------------ | ----------- | ------------------- | -------- | --------------------- |
| 1 | Cross-cutting | Rescale/level transition | Standard `ckks.Rescale/RescaleTo` | `MISSING_FAST_IMPLEMENTATION` | Full-Q dropped-prime rounding, INTT/NTT, all lower limbs | BLOCKER | Task 003 — Fast q0/q1 Rescale contract |
| 2 | EvalMod | NTT ciphertext arithmetic API | Standard generic polynomial evaluator | `FAST_EXISTS_BUT_NOT_WIRED` | Repeated full-Q Mul/MulRelin and QP relinearization | BLOCKER | Task 004 — NTT-resident Fast evaluator surface |
| 3 | CoeffsToSlots / SlotsToCoeffs | Factorized LinearTransform | Standard common LT | `FAST_EXISTS_BUT_NOT_WIRED` | QP decomposition, hoisted GadgetProducts, ModDown, many rotations | BLOCKER | Task 005 — Fast DFT adapter |
| 4 | ScaleDown / ModUp | q0/q1 modulus lifecycle | Full-Q Standard boundary | `BOUNDARY_ONLY` | Full INTT, q0->all-Q redistribution, full NTT, optional QP switch | BLOCKER | Task 006 — Fast ScaleDown/ModUp/Trace boundary |
| 5 | EvalMod | Mod1/Paterson-Stockmeyer wiring | Standard Mod1 and polynomial evaluators | `STANDARD_HIGH_COST` | Repeated Mul, Rescale, temporary powers | BLOCKER | Task 007 — Fast EvalMod integration |
| 6 | Pack/unpack | Sparse ciphertext packing | Full-Q ring operations | `MISSING_FAST_IMPLEMENTATION` | All active Q limbs and ciphertext copies | No for one ciphertext; yes for batching | Task 008 — q0/q1 pack/unpack |
| 7 | N1/N2 boundary | Ring-degree switch | Standard ApplyEvaluationKey | `FAST_EXISTS_BUT_NOT_WIRED` | GadgetProduct, QP ModDown, temporary ciphertext | Configuration-dependent | Task 008 — ring-degree adapter |
| 8 | Bootstrap orchestration | Evaluator selection and stage calls | Standard embedded evaluator | `MISSING_FAST_IMPLEMENTATION` | Hidden Standard fallback across all stages | BLOCKER | Task 009 — end-to-end Fast Bootstrap evaluator |
| 9 | CI support | Real/Complex switching | Standard DomainSwitcher | `UNKNOWN_REQUIRES_TARGETED_AUDIT` | Evaluation-key/QP ring-type switch | Configuration-dependent | Follow-up after Standard-ring Bootstrap |
| 10 | Hot-path cleanup | Scratch/domain reuse | Per-call temporaries/conversions | `MISSING_FAST_IMPLEMENTATION` | Allocations and repeated partial transforms | No | Task 010 — profile and remove measured debt |

## 16. Stage-A implementation roadmap

### Task 003 — Define and implement Fast q0/q1 Rescale

- Behavior: preserve q0/q1 authority, scale updates, and explicit level metadata without reading dormant limbs; define `Rescale` and `RescaleTo` contracts for NTT-resident Fast ciphertexts.
- Likely files: `schemes/ckks/fast/evaluator.go`, a new bounded Fast rescale file, and Fast-only tests.
- Reuse: q0/q1 subring operations and existing metadata conventions.
- Keep untouched: `schemes/ckks/evaluator.go` and `ring/scaling.go` Normal behavior.
- Minimum tests: one- and multi-level transitions, in-place/out-of-place, scale updates, dormant-limb non-access, decoded/authoritative reference behavior, rejection at insufficient level.
- Priority reason: DFT and EvalMod both repeatedly require it; no end-to-end path can preserve current circuit contracts without a level policy.

### Task 004 — Add an NTT-resident Fast evaluator surface

- Behavior: implement the Stage-A subset of `schemes.Evaluator` needed by DFT/Mod1: ciphertext Mul/MulRelin, truncation-backed Relinearize, scalar/plaintext Mul and MulThenAdd, Add/Sub, and Rescale dispatch, all on q0/q1.
- Likely files: `schemes/ckks/fast/evaluator.go` plus small operation-specific Fast files/tests.
- Reuse: FastAdd/FastSub, current Fast multiplication formulas, FastTruncateDegree2To1, Task 003 Rescale.
- Keep untouched: `ckks.Evaluator` and Standard RLWE GadgetProduct.
- Minimum tests: interface compile assertion, NTT/non-Montgomery Bootstrap representation, aliasing, scale/level metadata, polynomial power-basis operation sequences, dormant-limb poison tests.
- Priority reason: unlocks the generic polynomial call shape and avoids catastrophic per-Mul NTT/INTT wrapping.

### Task 005 — Wire a Fast DFT evaluator

- Behavior: execute each CoeffsToSlots/SlotsToCoeffs factor matrix with Fast LinearTransform/automorphism, q0/q1 additions/scalars, and Fast Rescale; support current factorization levels.
- Likely files: a Fast DFT adapter under `circuits/ckks/dft` or `circuits/ckks/bootstrapping`, plus narrowly scoped Fast LinearTransform contract changes.
- Reuse: Fast LinearTransform, FastAutomorphism/Rotate, FastAdd/FastSub, Tasks 003–004.
- Keep untouched: common Standard LinearTransform and its BSGS/QP implementation.
- Minimum tests: one factor, multi-factor level progression, CoeffsToSlots and SlotsToCoeffs semantic comparisons, dense/sparse formats, no QP/GadgetProduct path.
- Priority reason: removes the largest repeated QP rotation/ModDown region while reusing an existing primitive.

### Task 006 — Define and implement Fast ScaleDown/ModUp/Trace boundary

- Behavior: retain q0/q1 authority across the modulus-raising boundary, preserve required scale/metadata, and perform sparse Trace with Fast automorphisms; explicitly decide whether Fast ScaleDown stops at level 1.
- Likely files: Fast-only Bootstrap stage file(s) and targeted tests; reference helper use only in tests/debug.
- Reuse: q0/q1 reconstruction as oracle, Fast automorphism, Fast Add/Sub, Fast scalar operations.
- Keep untouched: Standard `bootstrapping.Evaluator.ModUp` and Standard dense/sparse key switching.
- Minimum tests: boundary scale/level contract, coefficient/NTT transitions, dormant-limb poison, sparse and dense slot settings, no full-Q redistribution in production path.
- Priority reason: unavoidable Bootstrap boundary with high structural cost; deferred until the evaluator level contract is explicit.

### Task 007 — Integrate Fast EvalMod

- Behavior: run current Mod1/Paterson-Stockmeyer structure through the Fast evaluator, including power basis, DoubleAngle, and optional inverse polynomial where configured.
- Likely files: a Fast Mod1/polynomial adapter or evaluator injection point plus integration tests.
- Reuse: Task 004 evaluator and Task 003 Rescale; do not redesign current Fast Mul.
- Keep untouched: Standard Mod1 and polynomial evaluators.
- Minimum tests: representative Mod1 parameter sets, expected level/scale progression, decoded semantic comparison, proof that Standard GadgetProduct is not reached.
- Priority reason: completes the multiplication-heavy middle of Bootstrap after its primitive dependencies exist.

### Task 008 — Wire q0/q1 packing and ring-degree conversion

- Behavior: pack/unpack only authoritative limbs and replace Standard N1/N2 ApplyEvaluationKey with existing Fast ring-degree conversion when parameters differ.
- Likely files: Fast Bootstrap packing adapter and, only if needed, bounded contract changes in `schemes/ckks/fast/ring_degree.go`.
- Reuse: FastN1ToN2/FastN2ToN1 and q0/q1 subring multiplication.
- Keep untouched: Standard packing and ApplyEvaluationKey paths.
- Minimum tests: one/many ciphertexts, odd batch count, N1=N2 and N1!=N2, NTT representation, metadata/log-slot restoration, dormant-limb poison.
- Priority reason: configuration-dependent coverage with existing core primitive; lower priority than single-ciphertext Bootstrap stages.

### Task 009 — Add the Fast Bootstrap orchestrator

- Behavior: select and compose Fast Stage-A implementations for the public Bootstrap flow while preserving the Normal evaluator and API-facing structure.
- Likely files: Fast-specific bootstrapping evaluator/orchestration files and integration tests.
- Reuse: Tasks 003–008 and existing Fast Bootstrap key support.
- Keep untouched: existing Normal `bootstrapping.Evaluator` behavior.
- Minimum tests: complete Standard-ring Bootstrap, dense/sparse slots, N1=N2 and N1!=N2, post-Bootstrap supported operation, decoded semantics, explicit assertion against hidden Standard/QP fallback.
- Priority reason: proves Stage-A coverage and enables end-to-end performance measurement.

### Task 010 — Profile and remove measured domain/allocation debt

- Behavior: benchmark Fast versus Normal by stage, then reuse scratch, retain NTT domain, and remove only measured hot allocations/conversions.
- Likely files: Fast evaluator/DFT/Bootstrap internals and benchmark tests.
- Reuse: completed Fast engine.
- Keep untouched: Normal implementation and Stage-B noise knobs.
- Minimum tests: unchanged correctness suite plus allocation/benchmark baselines.
- Priority reason: optimization follows complete coverage, with evidence from the runnable engine.

## 17. Recommended next task

Proceed with **Task 003 — Define and implement Fast q0/q1 Rescale**. It has the highest expected Bootstrap speed/coverage impact relative to implementation scope because it is shared by ScaleDown, both DFTs, the polynomial power basis, Paterson-Stockmeyer combination, and DoubleAngle. The task must begin by fixing the Fast mathematical contract for level consumption; it must not call Standard full-Q Rescale on dormant limbs.
