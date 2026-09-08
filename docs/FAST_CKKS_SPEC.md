# Fast-CKKS specification

This document is the single detailed source of truth for Fast-CKKS work on the `fast-ckks` branch. It records the current implementation, the permanent architectural boundaries, and the future experiments. Statements marked **Target design**, **Requires source audit**, or **Future experiment** are not claims that the current code already implements them.

## 1. Project objective and priority

The application keeps ordinary CKKS API usage while the backend supplies Normal or Fast behavior:

```text
Application
    |
    | unchanged CKKS API usage
    v
Lattigo
    |
    +-- upstream/normal behavior
    |
    +-- Fast-CKKS behavior on fast-ckks branch
```

The application must not contain Normal/Fast branches. Fast behavior belongs in the backend library.

The primary engineering goal of the current stage is speed. Fast-CKKS must make CKKS execution sufficiently fast for large workloads, including thousands of CNN inferences used to collect FHE numerical-error data. If Fast execution is not substantially faster than Normal Lattigo, the current project goal has not been met. The project is intentionally insecure; security is not a current optimization constraint.

Development priority:

1. Fast execution throughput.
2. Functional correctness of the simplified Fast execution.
3. Clean extension points for later experiments.
4. Noise fidelity evaluation and improvement later.

Current Fast execution is not required to reproduce the full Standard FHE noise distribution. Do not add expensive computation merely to make Fast noise look more realistic.

## 2. Non-goals and scope

Fast-CKKS is not intended to redesign CKKS or the application, modify CNN/model code, reduce the frontend-configured parameter chain, refactor all Lattigo schemes, create Fast BFV/BGV/TFHE, or optimize unrelated code.

Primary scope:

- `schemes/ckks` and its Fast backend;
- CKKS circuits required by the application;
- `circuits/ckks/bootstrapping`; and
- directly required `core/rlwe` execution paths.

Other schemes remain untouched unless a genuinely shared primitive is unavoidably required by a CKKS call path. Preserve the original Normal implementation; Fast code remains separate where practical so both paths support correctness comparison, performance comparison, and later noise calibration.

## 3. Core design principle

> Preserve structure; elide computation.

Preserve public APIs where practical, CKKS parameter objects, metadata, structural compatibility, Normal behavior, and enough representation structure for future experiments. In Fast hot paths, do not compute expensive unneeded residue values, perform full-RNS work merely because storage exists, perform security-only work merely because Standard Lattigo does, synchronize every stored limb per operation, or invoke Standard QP/key-switch paths as a hidden fallback.

Unused Fast-state fields or residue arrays may exist structurally while their values are unmaintained, stale, or dormant. The optimization target is to avoid computing, updating, transforming, or synchronizing RNS residue values that are not needed for the Fast numerical result.

## 4. Parameter and RNS representation

Fast and Normal execution use the same application-provided CKKS parameter configuration. A 20-level modulus chain remains a 20-level configuration; Fast must not reinterpret it as two levels.

Use the following terms precisely:

```text
q_i      = modulus prime at RNS index i
r_i[k]   = x[k] mod q_i, the stored residue for coefficient k in limb i
c0,c1... = ciphertext polynomial components
Level    = highest active logical Q index
```

When the current Level is at least 1, Fast execution may maintain only a small sufficient subset of RNS residue values—currently the first two limbs, r0 and r1—when those values are sufficient for the required Fast reconstruction or arithmetic shortcut. For a degree-1 ciphertext this can be pictured as:

```text
c0                         c1
 ├─ r0 actively maintained    ├─ r0 actively maintained
 ├─ r1 actively maintained    ├─ r1 actively maintained
 └─ r2...rL unmaintained      └─ r2...rL unmaintained
```

This is a computational shortcut, not a requirement that every Fast ciphertext physically retain Level >= 1. Fast follows Standard CKKS Level progression, including transitions to Level 0. At Level 0 only the residue limb for modulus q0, containing r0 values, naturally exists.

The configured modulus values `q0, q1, ..., qL` remain in the parameter object even when a hot path does not actively maintain all corresponding residue arrays:

```text
modulus parameter exists
!=
its corresponding residue values are actively maintained
```

Do not redefine historical q0/q1 Fast terminology as c0/c1. Historical “q0/q1-only” wording means an r0/r1 residue-maintenance shortcut, not deletion of modulus parameters or a floor on ciphertext Level.

## 5. Fast hot-path rules

In production Fast hot paths, avoid unless explicitly required at a deliberate boundary:

- CRT reconstruction;
- redistribution into unneeded `r2...rL` storage;
- full-Q NTT/INTT;
- full-Q arithmetic;
- QP arithmetic;
- GadgetProduct;
- Standard key switching or Standard relinearization;
- unnecessary basis extension or ModDown;
- unnecessary temporary allocations; and
- unnecessary coefficient-domain/NTT-domain round trips.

Reference/debug helpers such as r0/r1 reconstruction may remain available, but must not silently migrate into production hot paths. Dormant residue storage must not be read, synchronized, or treated as mathematically current merely because it exists. Follow Standard CKKS Level/Scale semantics when an operation requires a level transition.

### 5.1 Reconstruction representation contract

For coprime q0 and q1, `(r0, r1)` uniquely determines an element of `Z/(q0*q1)Z`. CRT uniqueness is not an architectural risk. The implementation-sensitive issue is normalization and interpretation.

Before ordinary coefficient reconstruction, the conceptual precondition is:

```text
coefficient domain
non-Montgomery representation
canonical residue r_i in [0, q_i)
```

- NTT coordinates require an appropriate INTT boundary before coefficient reconstruction.
- Montgomery values require conversion to normal representation.
- Lazy arithmetic may produce values in `[0, 2*q_i)`, `[0, 3*q_i)`, or another documented lazy range; reduce them to canonical residues before CRT interpretation.

Reference/debug code may reject nonconforming input. Production code may normalize efficiently when its caller contract permits. This requirement does not authorize routine NTT -> INTT -> NTT round trips: keep the Fast evaluator NTT-resident where possible and convert only when reconstruction is mathematically required.

### 5.2 ScaleDown semantics

For current Level `level`, define:

```text
Q_level = q0 * q1 * ... * q_level
Delta   = ciphertext Scale
rho     = Mod1 MessageRatio

currentMessageRatio = Q_level / Delta
```

Current Lattigo `ScaleDown` first removes unnecessary upper moduli while the source check is conceptually equivalent to:

```text
currentMessageRatio >= q_level * rho
```

which implies:

```text
(Q_level / q_level) / Delta >= rho
```

The loop uses `Resize`/DropLevel. It does not divide polynomial coefficients, change Scale, or introduce rounding; it only removes an unnecessary high-modulus representation. This is distinct from Rescale.

The three cases are:

1. **Already Level 0:** no level dropping and no `RescaleTo`; scalar adjustment may still occur.
2. **Free DropLevel reaches Level 0:** remove all unnecessary upper levels, apply scalar adjustment, and do not call `RescaleTo`.
3. **Free DropLevel stops at Level j > 0:** compute the current message ratio, multiply coefficients by an integer approximately `scaleUp = currentMessageRatio / targetMessageRatio`, update Scale by the same factor, then call `RescaleTo` to reach the target Bootstrap scale and eventually Level 0.

DropLevel is preferred first because it is cheaper and introduces no rescaling rounding error. If both coefficients and Scale are multiplied by `a`, the represented message is unchanged:

```text
x -> a*x
Delta -> a*Delta
(a*x)/(a*Delta) = x/Delta
```

The scalar adjustment aligns message-ratio/scale state with the Bootstrap target, whose final scale is approximately `q0 / MessageRatio` at Level 0.

### 5.3 Fast Rescale semantics

Fast Rescale preserves Standard CKKS mathematics. For one consumed highest logical modulus `q_L`:

```text
x'     = round(x / q_L)
Scale' = Scale / q_L
Level' = Level - 1
```

Current Lattigo may consume multiple moduli according to `LevelsConsumedPerRescaling()`. Fast changes how the result is computed, not what Rescale means.

Standard Lattigo maintains all active residues and uses the residue for the dropped modulus in `DivRoundByLastModulusManyNTT`. If Fast has not maintained that residue, it must not call Standard Rescale on stale data. When Level >= 1 and r0/r1 are the maintained shortcut residues, a Fast implementation may normalize them, reconstruct modulo `q0*q1`, choose the required centered representative, round-divide by logical divisor `q_L`, reduce the quotient into residues that remain actively maintained, and update Scale and Level exactly as Standard semantics require. A transition to Level 0 retains only r0; there is no artificial Level-1 floor.

The existing per-coefficient `big.Int` reconstruction is a reference oracle, not a production Rescale implementation. Production Stage A should investigate allocation-free fixed-width 64/128-bit reconstruction and rounded division, validate its supported modulus-size contract, avoid stale high residues, and leave Normal Rescale unchanged.

### 5.4 Target parameter profile

The supplied engineering profile, not a universal CKKS invariant, constructs 17 Q primes (`MaxLevel = 16`) in this order:

| Segment | Prime bit sizes |
|---|---|
| Q0 | 55 |
| SlotsToCoeffs Q | 39, 39, 39 |
| Circuit Q | 45 |
| EvalMod Q | eight 60-bit primes |
| CoeffsToSlots Q | four 56-bit primes |

Thus `q0 < 2^55`, `q1 < 2^39`, and `q0*q1 < 2^94`. For this profile, per-coefficient r0/r1 reconstruction does not require arbitrary-precision arithmetic and fits in a two-`uint64`/`math/bits` candidate representation. Future parameter sets must be validated rather than assumed to share this bound.

Production direction:

```text
Reference oracle:
    existing big.Int r0/r1 reconstruction + direct rounded division

Production candidate:
    allocation-free fixed-width r0/r1 reconstruction and rounded division
    using 64/128-bit arithmetic
```

## 6. Current implementation status

This branch does not start from Phase 0. The status below is derived from the current source and branch history through `a03cec77`.

| Area | Status | Current evidence and boundary |
|---|---|---|
| r0/r1 reconstruction and redistribution utilities | Reference/debug only | `schemes/ckks/fast/fast.go`; useful as a correctness oracle or explicit materialization boundary, not the production hot path |
| First-two-limb partial NTT/INTT | Implemented | `partial_ntt.go`; transforms only residue storage for q0 and q1 |
| Fast key generation | Implemented | `keys.go`; current key lifecycle uses zero-secret semantics and Fast layouts |
| First-two-limb evaluator arithmetic | Implemented | NTT-resident Add/Sub, Mul/MulRelin, scalar/plaintext Mul, fused MulThenAdd, truncation-backed Relinearize, and Rescale operate on maintained r0/r1 residues; full `schemes.Evaluator` conformance remains blocked by its unsafe full-Q/QP `rlwe.EvaluatorProvider` surface |
| Fast multiplication | Implemented | `FastMulQ01Authoritative` avoids CRT, bounds, and redistribution; reference multiplication remains available |
| Degree truncation | Implemented | `FastTruncateDegree2To1` under the current zero-secret mode |
| Fast automorphism/rotation | Implemented | First-two-limb Standard-ring path; no evaluation-key use in the helper |
| Ring-degree conversion | Implemented | r0/r1-maintaining N1/N2 conversion |
| Fast LinearTransform | Implemented | single-level first-two-limb path without evaluation-key, QP, or level transition |
| Fast Bootstrap key support | Implemented | `circuits/ckks/bootstrapping/fast_keys.go` generates Fast evaluation-key material |
| Fast Rescale | Implemented | `schemes/ckks/fast/rescale.go` provides fixed-width q0/q1 arithmetic with Standard Scale/Level semantics; targeted tests pass |
| ModUp/ModDown, KeySwitch, and full Bootstrap execution | Not implemented | These remain future bounded stages; do not imply that Fast key support is Fast Bootstrap execution |

The current implementation status is not a permanent architecture claim. In particular, current zero-secret behavior and current evaluation-key material are implementation modes that future experiments may extend.

## 7. Current implementation versus permanent architecture

### CURRENT IMPLEMENTATION

The current Fast key path forces the secret-key contribution to zero when generating Fast material. The Fast public-key target is `(e_pk, 0)`: the error component is retained and the second component is the zero polynomial. Current evaluation-key generation is approximately error-only / `(e, 0)`-style according to the implemented Fast key-generation path. Current degree-2 truncation relies on this zero-secret mode, and current Fast operations support a bounded subset of CKKS execution.

These are current implementation facts, not permanent scientific assumptions. Do not state that `s = 0` is permanently required or that all future key material must remain zero/error-only.

### PERMANENT FAST ARCHITECTURE

The durable boundaries are: preserve Normal behavior and public structure where practical; maintain only the minimum sufficient residue subset for each Fast operation; permit normal Level/Scale progression through Level 0; avoid hidden Standard/full-RNS fallback on stale residue storage; preserve the full CKKS parameter chain; and leave localized extension points for alternate secrets, key models, and equivalent-noise experiments.

## 8. Key policy and future experimental knobs

Fast mode intentionally removes or elides computational contributions of relevant key material in its current implementation. The exact semantics of future nonzero-secret evaluation keys require a dedicated source-level and mathematical audit.

Planned conceptual backend experiment knobs (not implemented here) include:

```text
SecretMode:
    Zero
    StandardSampled

PublicKeyErrorDistribution
PublicKeyErrorMagnitude
RandomSeed
EquivalentNoiseModel
```

`StandardSampled` should use a normal Lattigo secret polynomial drawn from the configured `Xs`, allowing paired Normal/Fast experiments to use the same sampled secret. Do not add application/frontend Fast flags; these knobs belong in a backend experiment harness or configuration layer.

Initially, later experiments should vary Public Key error only. Do not make evaluation-key error an experiment variable yet.

For future nonzero-secret experiments, do not assume that setting the `a` polynomial to zero means all secret-dependent terms may be removed, auxiliary P moduli may be removed, or Gadget/RNS decomposition has no relevance. The roles of P, secret-dependent gadget terms, key switching, relinearization, and Galois keys require explicit audit.

## 9. Two development stages

### Stage A — Fast Engine

Current priority: complete a runnable CKKS/Bootstrap execution path with maximum practical speed. Cleaner noise than Standard FHE is acceptable. Stage A should make large workloads such as thousands of CNN inferences practical. First eliminate asymptotically or structurally expensive Standard work, then complete Fast coverage; do not prematurely micro-optimize before the end-to-end path exists.

### Stage B — Noise Experiments

After the Fast Engine works end-to-end, use the experimental knobs and, if needed, theory-guided equivalent-noise injection to study FHE numerical error. Stage B extends Stage A rather than replacing it.

Removing security-related work may alter numerical-noise propagation in key switching, relinearization, rotation, ModDown, rescaling, and bootstrapping. A later operation-specific equivalent-noise model may inject the omitted effect efficiently. Do not add guessed Gaussian noise or implement this model during Stage A.

## 10. Operations and Bootstrap audit

Before implementing additional behavior, independently trace `Mul`, `MulRelin`, `Relinearize`, `Rotate`, `Rescale`, and `Bootstrap`. Record the CKKS entry point, downstream functions, ciphertext degree, components read/written, keys used, key-switch location, RNS/level behavior, implementation layer, and smallest modification point.

Bootstrap is the highest-risk integration stage. Trace the actual current Lattigo flow before implementation, including rotations, key switching, relinearization, multiplication, rescaling, modulus/level transitions, evaluation keys, and degree changes. Do not copy the complete Standard Bootstrap implementation merely to create a Fast version; prefer the smallest reusable operation boundary.

## 11. Validation terminology and methods

Use three domains precisely:

- **Message domain:** the original vector before CKKS encoding.
- **Plaintext domain:** the encoded CKKS polynomial/plaintext representation.
- **Ciphertext domain:** RLWE ciphertext polynomial components.

Do not call the pre-encoding vector a plaintext.

### Fast correctness validation

Use Normal Lattigo as a reference where mathematically appropriate. For individual Fast primitives, compare actively maintained residues with the corresponding Normal/reference operation when contracts are expected to match. For full Bootstrap, compare decoded semantic outputs rather than requiring raw ciphertext equality unless exact equality is justified.

### Noise experiment

Given `ct_before` and `ct_after = Bootstrap(ct_before)`, define an aligned difference:

```text
ct_error = ct_after - ct_before
```

Before subtraction, ensure compatible secret, ring degree, ciphertext degree, scale, level, and domain/NTT representation. Then use:

```text
message_error = Decode(Decrypt(ct_error))
```

for statistical analysis and later CNN robustness training. Retain a sanity check against:

```text
Decode(Decrypt(ct_after)) - Decode(Decrypt(ct_before))
```

within expected CKKS tolerance.

The target workload uses polynomial ring degree `N = 65536`; this is not 65,536 CKKS slots. Standard complex CKKS normally has at most `N/2` slots. Derive slot counts from the active parameter configuration rather than hard-coding them.

## 12. Performance and testing

Performance is a first-class acceptance criterion during development. Each major Fast primitive should eventually be measurable against its Normal equivalent. Do not defer all performance measurement until semantics are complete, but do not prematurely micro-optimize before the complete Fast execution path exists. After coverage is complete, profile allocation, domain-conversion, and scratch-buffer hotspots.

Testing should include:

1. Normal regression where Normal behavior remains available.
2. Fast invariant tests for actively maintained residues, Level/Scale semantics, and intentional zero-key states.
3. Independent tests for each modified primitive.
4. Bootstrap integration and continued supported use of its result.
5. Application compatibility without Fast-specific application changes.

## 13. Known Performance Questions

Non-binding questions for investigation after Fast coverage is complete:

- Can repeated temporary polynomial allocation use reusable buffers or pools?
- Can repeated first-two-residue-limb NTT/INTT conversions be reduced?
- Can Fast ciphertexts remain in NTT domain across larger operation sequences?
- Can automorphism scratch allocations be reused?
- Can LinearTransform and Bootstrap reuse scratch buffers?

These are not requests to optimize code in this documentation task.

## 14. Hard prohibitions

- [ ] Do not modify application/model architecture or add frontend Fast/Normal flags.
- [ ] Do not reduce the configured CKKS parameter chain to two levels.
- [ ] Do not modify BFV/BGV/TFHE merely for symmetry.
- [ ] Do not infer current Lattigo behavior from names alone.
- [ ] Do not confuse modulus primes `q_i`, stored residues `r_i[k]`, ciphertext components `c0/c1/...`, or ciphertext Level.
- [ ] Do not read or synchronize unmaintained residue storage in Fast hot paths without an explicit operation requirement.
- [ ] Do not impose an artificial Level-1 floor; follow Standard CKKS Level/Scale semantics through Level 0.
- [ ] Do not call Normal key-dependent/full-RNS paths with zero keys and label that implementation complete.
- [ ] Do not silently restore security semantics or add guessed noise during Stage A.
- [ ] Do not create repository-wide refactors.
- [ ] Do not start Bootstrap implementation before tracing its actual source path.

## 15. Definition of done

The overall project eventually satisfies:

```text
Same application source
        |
        +-- standard Lattigo -> Normal behavior
        |
        +-- Fast-CKKS fork/branch -> Fast behavior
```

The application remains unchanged. Fast-specific details remain in the backend. Stage A has a complete runnable fast path with measurable throughput advantage; Stage B can then add localized, theory-guided noise experiments without replacing the engine.
