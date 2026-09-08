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

Preserve public APIs where practical, CKKS parameter objects, metadata, structural compatibility, Normal behavior, and enough representation structure for future experiments. In Fast hot paths, do not compute expensive dormant values, perform full-RNS work merely because storage exists, perform security-only work merely because Standard Lattigo does, reconstruct or redistribute `q2...qL` per operation, or invoke Standard QP/key-switch paths as a hidden fallback.

Unused Fast-state fields or limbs may exist structurally while their values are dormant and non-authoritative.

## 4. Parameter and RNS representation

Fast and Normal execution use the same application-provided CKKS parameter configuration. A 20-level modulus chain remains a 20-level configuration; Fast must not reinterpret it as two levels.

The authoritative Fast representation is:

```text
q0, q1 = authoritative RNS limbs
q2 ... qL = dormant / non-authoritative during Fast hot-path execution
```

This is separate from ciphertext degree. Use `c0`, `c1`, `c2`, ... for ciphertext polynomial components and `q0`, `q1`, `q2`, ... for RNS modulus limbs. A degree-1 ciphertext can conceptually contain:

```text
c0                         c1
 ├─ q0 authoritative         ├─ q0 authoritative
 ├─ q1 authoritative         ├─ q1 authoritative
 └─ q2...qL dormant          └─ q2...qL dormant
```

Do not redefine historical q0/q1 Fast terminology as c0/c1. The full CKKS parameter chain remains configured even when only q0/q1 are authoritative.

## 5. Fast hot-path rules

In production Fast hot paths, avoid unless explicitly required at a deliberate boundary:

- CRT reconstruction;
- redistribution to `q2...qL`;
- full-Q NTT/INTT;
- full-Q arithmetic;
- QP arithmetic;
- GadgetProduct;
- Standard key switching or Standard relinearization;
- unnecessary basis extension or ModDown;
- unnecessary temporary allocations; and
- unnecessary coefficient-domain/NTT-domain round trips.

Reference/debug helpers such as q0/q1 reconstruction may remain available, but must not silently migrate into production hot paths. Dormant limbs must not be read, synchronized, or treated as mathematically current merely because storage exists.

## 6. Current implementation status

This branch does not start from Phase 0. The status below is derived from the current source and branch history through `a03cec77`.

| Area | Status | Current evidence and boundary |
|---|---|---|
| q0/q1 reconstruction and redistribution utilities | Reference/debug only | `schemes/ckks/fast/fast.go`; useful as a correctness oracle or explicit materialization boundary, not the authoritative hot path |
| q0/q1 partial NTT/INTT | Implemented | `partial_ntt.go`; transforms only q0 and q1 |
| Fast key generation | Implemented | `keys.go`; current key lifecycle uses zero-secret semantics and Fast layouts |
| q0/q1 evaluator arithmetic | Implemented | Fast evaluator and arithmetic paths operate on authoritative limbs |
| Fast multiplication | Implemented | `FastMulQ01Authoritative` avoids CRT, bounds, and redistribution; reference multiplication remains available |
| Degree truncation | Implemented | `FastTruncateDegree2To1` under the current zero-secret mode |
| Fast automorphism/rotation | Implemented | q0/q1-only Standard-ring path; no evaluation-key use in the helper |
| Ring-degree conversion | Implemented | q0/q1-authoritative N1/N2 conversion |
| Fast LinearTransform | Implemented | single-level q0/q1-authoritative path without evaluation-key, QP, or level transition |
| Fast Bootstrap key support | Implemented | `circuits/ckks/bootstrapping/fast_keys.go` generates Fast evaluation-key material |
| Fast Rescale, ModUp/ModDown, KeySwitch, and full Bootstrap execution | Not implemented | These remain future bounded stages; do not imply that Fast key support is Fast Bootstrap execution |

The current implementation status is not a permanent architecture claim. In particular, current zero-secret behavior and current evaluation-key material are implementation modes that future experiments may extend.

## 7. Current implementation versus permanent architecture

### CURRENT IMPLEMENTATION

The current Fast key path forces the secret-key contribution to zero when generating Fast material. The Fast public-key target is `(e_pk, 0)`: the error component is retained and the second component is the zero polynomial. Current evaluation-key generation is approximately error-only / `(e, 0)`-style according to the implemented Fast key-generation path. Current degree-2 truncation relies on this zero-secret mode, and current Fast operations support a bounded subset of CKKS execution.

These are current implementation facts, not permanent scientific assumptions. Do not state that `s = 0` is permanently required or that all future key material must remain zero/error-only.

### PERMANENT FAST ARCHITECTURE

The durable boundaries are: preserve Normal behavior and public structure where practical; keep q0/q1 authoritative; keep q2...qL dormant in Fast hot paths; avoid hidden Standard/full-RNS fallback; preserve CKKS parameters; and leave localized extension points for alternate secrets, key models, and equivalent-noise experiments.

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

Use Normal Lattigo as a reference where mathematically appropriate. For individual Fast primitives, compare authoritative q0/q1 results with the corresponding Normal/reference operation when contracts are expected to match. For full Bootstrap, compare decoded semantic outputs rather than requiring raw ciphertext equality unless exact equality is justified.

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
2. Fast invariant tests for authoritative limbs and intentional zero-key states.
3. Independent tests for each modified primitive.
4. Bootstrap integration and continued supported use of its result.
5. Application compatibility without Fast-specific application changes.

## 13. Known Performance Questions

Non-binding questions for investigation after Fast coverage is complete:

- Can repeated temporary polynomial allocation use reusable buffers or pools?
- Can repeated q0/q1 NTT/INTT conversions be reduced?
- Can Fast ciphertexts remain in NTT domain across larger operation sequences?
- Can automorphism scratch allocations be reused?
- Can LinearTransform and Bootstrap reuse scratch buffers?

These are not requests to optimize code in this documentation task.

## 14. Hard prohibitions

- [ ] Do not modify application/model architecture or add frontend Fast/Normal flags.
- [ ] Do not reduce the configured CKKS parameter chain to two levels.
- [ ] Do not modify BFV/BGV/TFHE merely for symmetry.
- [ ] Do not infer current Lattigo behavior from names alone.
- [ ] Do not treat q0/q1 as ciphertext components; use c0/c1 for those.
- [ ] Do not read or synchronize dormant q2...qL in Fast hot paths without an explicit boundary design.
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
