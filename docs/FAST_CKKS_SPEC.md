# Fast-CKKS specification

This document is the detailed source of truth for Fast-CKKS work on the `fast-ckks` branch. It describes the target design and the investigations required before implementation. Statements labeled **Requires source audit** are not confirmed behavior of the current Lattigo source.

## A. Project objective

The application keeps its ordinary CKKS API usage while the backend supplies normal or Fast behavior:

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

The application must not contain Normal/Fast mode branches. Fast behavior belongs in the backend library.

## B. Non-goals

Fast-CKKS is not intended to:

- provide cryptographic security in Fast mode;
- redesign CKKS or the application;
- modify CNN/model code;
- reduce the frontend-configured CKKS parameter chain;
- refactor all Lattigo schemes;
- create Fast versions of BFV, BGV, or TFHE; or
- optimize unrelated Lattigo code.

## C. Scheme scope

Primary scope:

- `schemes/ckks`;
- CKKS circuits required by the application;
- `circuits/ckks/bootstrapping`; and
- directly required `core/rlwe` execution paths.

Other scheme implementations remain untouched unless a genuinely shared primitive is unavoidably required. Any shared-code modification must be justified by an actual CKKS source-level dependency.

## D. Parameter invariant

Fast and Normal execution use the same application-provided CKKS parameter configuration. If the application configures a modulus chain corresponding to 20 levels, Fast mode receives that same configuration. Fast mode must not reinterpret this as `20 levels -> 2 levels`.

Parameter count/modulus-chain level and ciphertext degree are different concepts.

## E. Ciphertext terminology and representation

Lattigo RLWE ciphertexts have a ciphertext degree and polynomial components. For ordinary degree-1 ciphertexts, refer to the relevant components as `c0` and `c1`.

**Requires source audit:** verify the exact current Lattigo data structures before changing code. The historical Fast representation was described as “q0/q1-only”; implementation documentation must translate that intent into precise Lattigo terminology. Do not use `q0/q1` ambiguously when the source uses `q_i` for RNS moduli.

## F. Key policy

**Target design:** Fast mode removes the computational contribution of relevant key material. At minimum, investigate:

| Key material | Required investigation |
|---|---|
| Public key | Representation, lifecycle, and consumers |
| Relinearization key | Key-switching and degree-reduction paths |
| Galois keys | Rotation and permutation paths |
| Bootstrap/evaluation keys | Actual CKKS bootstrap call graph and stages |

**Requires source audit:** determine the exact zeroing location and lifecycle from current source. Do not assume all key types share a representation. Do not add validation that rejects the intentional Fast zero-key state.

Do not implement key zeroing as part of this documentation task.

## G. Execution objective

**Target design:** preserve the application’s CKKS API and parameter configuration while replacing expensive key-dependent backend execution with the intended simplified Fast execution. The application should continue using ordinary operations, including:

- Encrypt and Encode as applicable;
- Mul and MulRelin;
- Relinearize;
- Rotate;
- Rescale; and
- Bootstrap.

The application should not need APIs such as `FastMul`, `FastRotate`, or `FastBootstrap` unless source analysis proves an internal-only distinction necessary.

## H. Operations requiring source audit

Independently trace `Mul`, `MulRelin`, `Relinearize`, `Rotate`, `Rescale`, and `Bootstrap`. For every operation record:

- CKKS public entry point;
- downstream functions;
- ciphertext degree before and after;
- components read and written;
- evaluation keys used;
- key-switching location;
- RNS/level behavior;
- whether the implementation is CKKS-specific or shared RLWE; and
- the minimum possible modification point.

No implementation decision is valid purely from API naming.

## I. Bootstrap-specific investigation

Bootstrap is the highest-risk integration stage and should be implemented only after relevant primitive operations are understood. **Requires source audit:** trace the actual current Lattigo flow, including:

- major bootstrap stages;
- rotations and key switching;
- relinearization and multiplications;
- rescaling;
- modulus/level transitions;
- evaluation keys; and
- ciphertext-degree changes.

Do not copy the complete Bootstrap implementation merely to create a Fast version. Prefer the smallest reusable operation boundary.

## J. Proposed development phases

This is a planning framework, not a frozen implementation design.

### Phase 0 — Source audit

Read-only investigation producing:

- CKKS call graph;
- key-dependency map;
- ciphertext-component map;
- bootstrap map; and
- a `MUST MODIFY` / `MAY MODIFY` / `MUST NOT MODIFY` file list.

### Phase 1 — Fast key lifecycle

Implement and test only the verified key-zeroing lifecycle.

### Phase 2 — Primitive execution

Handle operations individually from audit results, likely including Mul/MulRelin, Relinearize, Rotate, and Rescale. Do not assume every operation needs a Fast implementation. Reuse normal code when analysis proves it is independent of removed key semantics and compatible with Fast invariants.

### Phase 3 — Bootstrap integration

Integrate Fast semantics into Bootstrap only after primitive operations are verified.

### Phase 4 — Integration validation

Run the unchanged application against the Fast Lattigo backend.

## K. Testing strategy

Maintain these categories:

1. **Normal regression:** ordinary CKKS behavior remains unchanged where Normal behavior remains available.
2. **Fast invariant tests:** intentional Fast states, including zeroed relevant key material, are verified.
3. **Operation tests:** each modified primitive is tested independently.
4. **Bootstrap integration:** Fast Bootstrap behavior is tested, including continued use of its resulting ciphertext in supported operations.
5. **Application compatibility:** the external application requires no Fast-specific source changes.

## L. Performance objective

Measure performance only after functional Fast semantics are established. Compare Normal operation time, Fast operation time, and speedup. Do not mix profiling implementation into the core Fast execution path unnecessarily.

## M. Hard prohibitions

- [ ] Do not modify application/model architecture.
- [ ] Do not introduce frontend Fast/Normal branches.
- [ ] Do not reduce CKKS parameter levels to two.
- [ ] Do not modify BFV/BGV/TFHE merely for symmetry.
- [ ] Do not infer Lattigo behavior from names alone.
- [ ] Do not call Normal key-dependent paths with zero keys and label that implementation complete.
- [ ] Do not silently restore security semantics.
- [ ] Do not create repository-wide refactors.
- [ ] Do not start Bootstrap implementation before tracing its actual source path.

## N. Definition of done for the overall project

Eventually, the same application source must work with both backends:

```text
Same application source
        |
        +-- standard Lattigo -> Normal behavior
        |
        +-- Fast-CKKS fork/branch -> Fast behavior
```

The application remains unchanged, and Fast-specific implementation details remain in the Lattigo backend.
