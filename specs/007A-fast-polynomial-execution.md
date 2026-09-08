# Task 007A — Fast Polynomial Execution Boundary

## Goal

Add a bounded CKKS-specific Fast polynomial evaluation path that can execute the polynomial stage required by Bootstrap EvalMod without exposing the Standard `schemes.Evaluator` / full-Q / QP / key-switching surface.

This task is the execution foundation for Task 007B Fast EvalMod.

Starting point:

```text
Fast ModUp + Trace        implemented
Fast CoeffsToSlots        implemented
Fast q0/q1 arithmetic     implemented
Fast Rescale              implemented
Fast EvalMod              not implemented
```

Base commit:

```text
e85a2cc5005ee8cded360b992707e29181be8eeb
feat(ckks/fast): implement Fast Trace and ModUp
```

Task 007A must produce:

```text
Fast CKKS ciphertext
    -> Fast-aware polynomial execution
    -> same polynomial mathematics / level-scale plan as current CKKS polynomial circuit
    -> q0/q1 authoritative result
    -> no Standard evaluator fallback
```

Do **not** wire Bootstrap EvalMod in this task. Do not implement the Mod1 cosine/double-angle wrapper yet.

---

## 1. Source-audit facts that define this task

Read the current source before editing, especially:

```text
circuits/ckks/mod1/mod1_evaluator.go
circuits/ckks/polynomial/polynomial_evaluator.go
circuits/common/polynomial/polynomial_evaluator.go
circuits/common/polynomial/power_basis.go
schemes/schemes.go
schemes/ckks/fast/evaluator_ntt.go
circuits/ckks/dft/fast.go
```

The current call path is conceptually:

```text
Mod1 Evaluate
    -> ckks/polynomial.Evaluator
        -> common/polynomial.Evaluator
            -> NewPowerBasis(ct, basis)
            -> GenPower(...)
            -> Paterson-Stockmeyer decomposition/execution
```

Important current facts:

1. `common/polynomial.NewPowerBasis` begins with `ct.CopyNew()`.
2. Generic polynomial execution is typed against `schemes.Evaluator`.
3. `schemes.Evaluator` embeds `rlwe.EvaluatorProvider`, which exposes Standard decomposition/QP/Gadget/key-switching operations.
4. The Fast evaluator deliberately does **not** implement `schemes.Evaluator` because doing so would create an unsafe hidden Standard fallback surface.
5. The Fast evaluator already provides the arithmetic primitives needed by polynomial execution: Add/Sub, Mul/MulRelin, MulThenAdd, Relinearize, Rescale, scalar arithmetic, and the current Level/Scale semantics.

Do not solve Task 007A by making `schemes/ckks/fast.Evaluator` implement the full `schemes.Evaluator` interface.

---

## 2. Preserve polynomial mathematics; replace only the execution engine

Task 007A is **not** a new polynomial approximation design.

Preserve the current CKKS polynomial planning/math where practical:

```text
polynomial coefficients
basis type
Chebyshev recurrence
Paterson-Stockmeyer decomposition
level schedule
scale schedule
```

The intended architecture is analogous to the existing Fast DFT adapter:

```text
Standard CKKS polynomial path
    generic Standard evaluator

Fast CKKS polynomial path
    bounded Fast adapter
    -> explicit *fastckks.Evaluator
    -> q0/q1 arithmetic only
```

Prefer adding a CKKS-specific Fast path under:

```text
circuits/ckks/polynomial
```

for example a `FastEvaluator`, rather than modifying scheme-agnostic common polynomial infrastructure.

Small reusable helpers may be added to `schemes/ckks/fast` only when they are genuinely arithmetic/storage primitives, not Mod1-specific logic.

---

## 3. Required public/internal Fast surface

Provide a bounded Fast polynomial API sufficient for the current Bootstrap caller. A shape such as the following is preferred:

```go
type FastEvaluator struct {
    // CKKS params
    // explicit *fastckks.Evaluator
    // reusable Fast polynomial workspace
}

func NewFastEvaluator(params ckks.Parameters, eval *fastckks.Evaluator) *FastEvaluator

func (eval *FastEvaluator) Evaluate(
    ct *rlwe.Ciphertext,
    p bignum.Polynomial,
    targetScale rlwe.Scale,
) (*rlwe.Ciphertext, error)
```

Equivalent bounded APIs are acceptable if they keep the same architecture.

Required scope for 007A:

```text
single polynomial
Chebyshev basis
Standard ring
NTT + Montgomery input
Degree-1 input
Level >= 1
q0/q1 authoritative
```

This is sufficient for the current default Bootstrap `CosDiscrete` Mod1 polynomial.

Optional only if essentially free:

```text
Monomial basis
```

Do not add PolynomialVector/mapping/masking support merely for generic completeness unless current EvalMod source actually requires it.

---

## 4. Reuse the existing Paterson-Stockmeyer planning logic where safe

Do not re-invent the polynomial decomposition mathematics if the existing exported data structures/planning functions can be reused safely.

The Fast path may reuse current common polynomial types and the current CKKS simulation/planning logic to derive:

```text
optimal split
Paterson-Stockmeyer sub-polynomials
required target Levels
required target Scales
```

Because a new Fast evaluator is added inside `circuits/ckks/polynomial`, it may reuse package-local CKKS scale simulation helpers where appropriate.

The unsafe boundary is **execution through `schemes.Evaluator`**, not the pure polynomial planning data.

Do not modify Standard polynomial results or Standard evaluator behavior.

---

## 5. Fast PowerBasis generation

Do not call the existing:

```go
common/polynomial.NewPowerBasis(ct, basis)
```

on the Fast production path if doing so triggers a full `ct.CopyNew()` of dormant high limbs.

Implement a Fast-aware PowerBasis/workspace whose authoritative numerical state is q0/q1 only.

For Chebyshev basis preserve the same recurrence used by current source:

```text
T_n = 2*T_a*T_b - T_|a-b|
where n = a + b
```

Preserve the same `SplitDegree` behavior / power-generation dependency structure unless source evidence requires otherwise.

Power generation must use the existing Fast arithmetic semantics:

```text
Fast Mul / MulRelin
Fast Add / Sub
Fast Relinearize where required
Fast Rescale
```

Never call Standard `MulRelin`, Standard `Relinearize`, Standard `Rescale`, GadgetProduct, QP arithmetic, or key switching.

---

## 6. Storage and allocation contract

The current `ring.Poly` / `rlwe.Ciphertext` representation ties logical Level to the number of coefficient rows. Therefore Task 007A does **not** require a global sparse-ciphertext redesign and must not pretend that a Level-12 ciphertext can structurally contain only two rows.

However, Fast must preserve this distinction:

```text
row exists structurally
!=
row is numerically authoritative
!=
row must be copied/computed every operation
```

Required behavior:

- q0/q1 are the only authoritative residues read and numerically updated.
- q2...qL must never influence q0/q1 output.
- Do not copy q2...qL merely because a source ciphertext has those rows.
- Do not perform full-Q arithmetic on q2...qL.
- Do not deliberately synchronize dormant rows.
- Metadata and logical Level/Scale semantics remain valid.

A structural allocation of high row storage caused by the current ciphertext representation is acceptable when unavoidable, but repeated warm-path allocation/copy of that storage should be avoided where practical.

---

## 7. Reusable workspace instead of `CopyNew`/`MulNew` churn

This project targets thousands of repeated evaluations. Do not build the Fast polynomial path around repeated full ciphertext `CopyNew()` calls.

The Fast polynomial evaluator should own or otherwise reuse workspace for repeatedly needed ciphertext values.

At minimum:

1. The `X^1`/`T_1` input copy must copy only q0/q1 + metadata into Fast-controlled storage.
2. Required power buffers should be lazily allocated and reusable across repeated evaluations of compatible parameters/degrees.
3. Avoid `MulNew` / `MulRelinNew` in the hot execution loop when an explicit reusable output buffer can be supplied instead.
4. Baby-step / giant-step temporary ciphertext storage should be reused where practical instead of forcing a fresh full-Level allocation for every sub-polynomial on every call.

It is acceptable for the first/cold evaluation to populate reusable workspace. Benchmark the warm steady state separately.

Do not create global mutable workspace. Evaluator-owned workspace is preferred and, like the existing Fast evaluator, may be documented as single-stream / not thread-safe.

---

## 8. Aliasing and metadata

The Fast polynomial evaluator may return evaluator-owned output or a newly allocated public result, but the ownership contract must be explicit and tests must prevent accidental overwrite on the next call.

Prefer a public API whose returned ciphertext remains valid after the call.

Preserve/produce correct:

```text
Degree
Level
Scale
IsNTT
IsMontgomery
IsBatched / LogDimensions metadata as applicable
```

The output Level and Scale must follow the same CKKS polynomial planner/simulator semantics as the current Standard polynomial circuit for the same input level, polynomial, and target scale.

Do not require raw q2...qL equality with Standard.

---

## 9. Validation contract

Reject invalid input rather than silently falling back.

Validate at least:

```text
nil evaluator / input / metadata
wrong ring type
input Degree != 1
Level < 1
non-NTT input
non-Montgomery input
insufficient levels for polynomial depth
unsupported polynomial basis
invalid/empty polynomial
```

If current source supports a zero-degree polynomial through a special path and it is cheap to preserve, support it; otherwise document the bounded requirement used by EvalMod.

Never make an invalid Fast case succeed by calling the Standard evaluator.

---

## 10. Correctness oracle strategy

Raw Fast ciphertext equality against Standard full-Q polynomial evaluation is not the only valid oracle because Fast deliberately leaves q2...qL dormant and uses its own maintained-residue execution semantics.

Use multiple focused oracles.

### A. Low-degree exact execution oracle

For small Chebyshev polynomials, build a test-only direct evaluation from already-tested Fast primitives and compare exact q0/q1 residues, Level, Scale, and representation.

The direct reference should be deliberately simple enough to catch errors in Fast PowerBasis / Paterson-Stockmeyer scheduling.

### B. Cleartext numerical oracle

Construct a deterministic input whose decoded numerical values are known under the current Stage-A zero-secret model, evaluate a polynomial, and compare the decoded q0/q1-authoritative result against clear evaluation of the same polynomial within a precision tolerance appropriate for CKKS rescaling.

Include a representative degree-30 Chebyshev polynomial, preferably the actual/default Mod1 polynomial generated from current bootstrapping parameters, without yet performing the Mod1 double-angle wrapper.

### C. Planner/metadata oracle

Verify output Level and Scale against the current CKKS polynomial simulation/planning result for the same degree/target scale.

Do not compare dormant high limbs.

---

## 11. Poisoned dormant-residue tests

Poison q2...qL before Fast polynomial evaluation.

Required result:

```text
q0/q1 result unchanged versus the same input with different high-limb poison
Level/Scale/metadata unchanged except for intended polynomial progression
q2...qL do not influence the answer
```

Where workspace outputs contain high rows, they may remain stale/uninitialized/dormant; do not require synchronization for raw ciphertext equality.

A test must fail if the Fast implementation accidentally routes an arithmetic operation through full-Q Standard code that consumes poisoned residues.

---

## 12. Required functional tests

Add focused tests covering at least:

### A. Small Chebyshev polynomial

Use a low degree where a simple direct recurrence/reference is easy to inspect.

Verify exact q0/q1 and Level/Scale.

### B. Power generation

Exercise powers requiring recursive split and Chebyshev recurrence, including at least one non-power-of-two degree.

Verify exact q0/q1 against direct Fast reference.

### C. Paterson-Stockmeyer path

Use a degree large enough to execute both baby and giant steps.

Verify against the direct/numerical oracle.

### D. Representative degree 30

Evaluate a degree-30 Chebyshev polynomial representative of Bootstrap Mod1.

Verify numerical output and expected level consumption.

### E. Poisoned high residues

Use two otherwise-identical inputs with different q2...qL poison and verify identical q0/q1 outputs.

### F. Repeated evaluation / workspace reuse

Run the same evaluator repeatedly with different q0/q1 inputs of the same shape.

Verify:

```text
no stale prior-input contamination
correct outputs
workspace reused
returned/public output ownership remains correct
```

### G. Validation

Exercise the invalid cases listed above.

---

## 13. Performance tests

Benchmark representative Fast polynomial evaluation at:

```text
LogN = 13
LogN = 16
```

Use a degree-30 Chebyshev polynomial and a Level sufficient for its planned depth.

Report separately:

```text
cold/first call (optional but useful)
warm steady-state
```

Warm steady-state must report:

```text
ns/op
B/op
allocs/op
```

The acceptance target is architectural rather than a hard speed ratio in 007A:

```text
no full-Q arithmetic
no QP/Gadget/key-switch work
no per-call full ciphertext CopyNew of X
no per-coefficient big.Int
q2...qL poison cannot affect q0/q1
warm allocations do not scale from repeatedly copying every dormant limb for every PowerBasis value
```

If structural `ckks.NewCiphertext` allocation remains a measurable dominant cost because of the current Level/storage representation, report it explicitly rather than hiding it. Do not redesign `ring.Poly` in this task.

A Standard full-Q polynomial benchmark may be reported for context, but do not claim a Fast/Standard ratio unless the compared paths evaluate the same polynomial with comparable setup/warm-state semantics.

---

## 14. Preserve existing tests

Run at least:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/polynomial
go test ./circuits/ckks/dft
go test ./circuits/ckks/bootstrapping
```

Do not regress:

```text
Fast Add/Sub/Mul/MulRelin/MulThenAdd
Fast Relinearize
Fast Rescale
Fast DFT
Fast ScaleDown
Fast ModUp/Trace
```

---

## 15. Out of scope

Do not implement in 007A:

```text
Mod1 / EvalMod wrapper
DoubleAngle loop
ArcSine / Mod1InvPoly integration
Bootstrap EvalMod wiring
full Bootstrap orchestration
PolynomialVector / slot mapping unless proven required by current EvalMod
DenseToSparse / SparseToDense
QP arithmetic
GadgetProduct
Standard key switching
full schemes.Evaluator conformance
new sparse ciphertext representation
ring.Poly / rlwe.Element Level redesign
per-coefficient arbitrary-precision CRT
noise-fidelity injection
```

Do not modify BFV/BGV behavior or Standard CKKS polynomial behavior.

---

## 16. Documentation

Update `docs/FAST_CKKS_SPEC.md` only after the implementation exists.

Add a status entry conceptually stating:

```text
Fast polynomial execution: implemented for the bounded single-polynomial CKKS path required by EvalMod
Fast EvalMod: not yet implemented
```

Also correct the existing ambiguous status wording if still present:

```text
Fast ModDown, KeySwitch, EvalMod, full Bootstrap execution | Not implemented
```

rather than implying the already-implemented Stage-A Fast ModUp boundary is still absent.

Do not claim EvalMod completion in 007A.

---

## Preferred files

Prefer a narrow change set around:

```text
circuits/ckks/polynomial/fast.go
circuits/ckks/polynomial/fast_test.go
circuits/ckks/polynomial/fast_bench_test.go   (or existing test file)
docs/FAST_CKKS_SPEC.md
```

Small changes to:

```text
schemes/ckks/fast
```

are allowed only if a missing narrow arithmetic/storage helper is required.

Avoid modifying:

```text
circuits/common/polynomial
schemes/schemes.go
core/rlwe
```

unless source evidence proves a tiny shared change is unavoidable. If so, explain why before broadening the change.

---

## Completion

Commit with:

```text
feat(ckks/fast): add Fast polynomial evaluator
```

Push to:

```text
origin/fast-ckks
```

Final report should include:

1. Fast polynomial API/location;
2. exact supported polynomial basis/scope;
3. how Paterson-Stockmeyer planning is reused;
4. Fast PowerBasis representation;
5. Chebyshev recurrence implementation;
6. confirmation that full `schemes.Evaluator` conformance was not added;
7. confirmation no Standard/QP/Gadget/key-switch fallback exists;
8. q0/q1-only copy/arithmetic strategy;
9. workspace reuse strategy;
10. high-limb poison result;
11. small-degree exact-reference result;
12. degree-30 numerical result;
13. output Level/Scale result;
14. repeated-evaluation/workspace-reuse result;
15. LogN13 warm benchmark;
16. LogN16 warm benchmark;
17. B/op and allocs/op;
18. any remaining structural allocation limitation;
19. tests run;
20. files changed;
21. commit hash;
22. push result.
