# Task 006B — Implement Fast ModUp Basis Raise

## Goal

Implement the Fast Bootstrap ModUp **basis-raise and scale-alignment portion only**.

Input is the valid output of Fast ScaleDown:

```text
Level = 0
NTT domain
only r0 exists logically
```

Output must be structurally compatible with the current Lattigo Bootstrap/DFT code:

```text
Level = BootstrappingParameters.MaxLevel()
NTT domain
r0 and r1 numerically maintained
r2...rL structurally present but dormant/stale
```

Do **not** implement Trace in this task. Task 006C will add Fast Trace and the final public ModUp behavior.

Base commit:

```text
132fc02a9144e671c8a359a804755d34513c2d77
perf(ckks/fast): specialize level-one Rescale
```

## Source behavior to preserve

Read the current Standard implementation in:

```text
circuits/ckks/bootstrapping/evaluator.go
```

For the Stage-A target configuration with no Dense↔Sparse evaluation keys, Standard ModUp does conceptually:

```text
Level-0 NTT ciphertext
-> INTT modulo q0
-> Resize to MaxLevel
-> centered lift from q0 into q1...qL
-> NTT
-> optional integer message scale-up
-> Trace
```

This task implements everything above **except Trace** and computes only the Fast-maintained q0/q1 numerical state.

## Important representation decision

Do not introduce a new sparse ciphertext type and do not modify `ring.Poly`, `rlwe.Element`, or the meaning of `Level()`.

Current Lattigo defines polynomial/ciphertext Level from the number of coefficient rows. Therefore a true:

```text
logical Level = MaxLevel
physical rows = only q0/q1
```

representation would require a broad core redesign and is out of scope.

For Task 006B:

```text
logical Level = MaxLevel
structural rows q0...qL exist
only r0/r1 are maintained
r2...rL are dormant/stale
```

This follows the existing Fast architecture: preserve structure, elide computation.

## Preserve backing storage when possible

Fast ScaleDown shrinks `ring.Poly.Coeffs` by slicing. The backing capacity and previously allocated high-level rows can still exist.

When raising Level back to MaxLevel, do not blindly call the existing `Poly.Resize(MaxLevel)` if reusable backing rows are still available.

Implement a narrow Fast helper that:

1. checks whether each component has capacity for `MaxLevel+1` coefficient rows;
2. if the hidden rows still have valid `[]uint64` backing arrays of length `N`, restores them by reslicing;
3. allocates only missing rows when the input genuinely lacks reusable backing storage.

Do not change global `ring.Poly.Resize` behavior.

The restored `r2...rL` contents are stale by definition and must not be read as current values.

## Input contract

Primary Stage-A input:

```text
ring.Standard
Degree = 1
Level = 0
NTT domain
ordinary / non-Montgomery representation
```

Reject Montgomery input for this task rather than adding conversion overhead.

Do not broaden to ConjugateInvariant, ring-degree switching, or other Bootstrap modes.

## Fast q0 -> q1 lift

At Level 0, only q0 is available.

For each ciphertext component:

1. INTT only q0 into coefficient representation.
2. Interpret each q0 residue using the same centered convention as the current Standard ModUp source.
3. Preserve q0 coefficient residue.
4. Compute only the corresponding q1 residue.
5. Do not compute q2...qL.
6. Forward NTT only q0 and q1.

Match the current Standard branch exactly for the q0-centered sign convention. Do not substitute a different centered-interval convention merely because it looks more symmetric.

Do not use per-coefficient `big.Int`.

Use fixed-width modular arithmetic / existing Ring helpers.

## Dormant high residues

After structural Level restoration:

```text
r2...rL
```

must not influence q0/q1 output.

Do not:

```text
full-Q INTT
full-Q NTT
full-Q centered redistribution
QP basis extension
GadgetProduct
ApplyEvaluationKey
Standard key switching
```

If high rows already contain poison/stale data, leave them dormant. Do not synchronize them merely to make storage look clean.

## Scale alignment

After q0/q1 NTT lift, preserve the same Standard scale-alignment rule:

```go
scale := (Mod1Parameters.ScalingFactor().Float64() / Mod1Parameters.MessageRatio()) / ct.Scale.Float64()
```

If:

```text
scale > 1
```

then:

1. `scalar = round(scale)` as in Standard;
2. multiply only maintained q0/q1 residues by this integer scalar;
3. update metadata scale using the same logical Standard value:

```text
ct.Scale = ct.Scale * scale
```

Do not accidentally set metadata Scale using `scalar` if Standard uses the floating `scale` value.

Reuse `MulIntegerMaintained` if appropriate, but do not broaden its contract unnecessarily.

If `scale <= 1`, leave ciphertext residues and Scale unchanged.

## API boundary

Extend the existing Bootstrap `FastEvaluator` in:

```text
circuits/ckks/bootstrapping
```

Prefer an internal/private method for this task, conceptually:

```go
modUpBasis(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error)
```

or an equivalently narrow name.

Do not expose a final public `ModUp` that claims full Standard ModUp semantics yet, because Standard ModUp includes Trace and Task 006C has not implemented Fast Trace.

Tests in the same package may exercise the private boundary directly.

## Standard/reference oracle

Because Standard public `ModUp` immediately performs Trace, create a small test/reference helper that mirrors the current **no-key Standard ModUp basis-raise + scale step** but stops before Trace.

The reference may use full-Q work because it is test/reference code.

Compare Fast against this reference for:

```text
q0 NTT residue
q1 NTT residue
Level
Scale
metadata/domain flags
```

Do not require equality of dormant q2...qL rows.

## Tests

Add focused tests for at least:

### A. Positive and negative centered coefficients

Use nontrivial values around the q0 centered boundary and compare q0/q1 exactly to the Standard-style reference.

### B. Scale-up path

Choose Scale so Standard's `scale > 1` branch executes.

Verify:

```text
integer residue multiplication
Scale metadata
q0/q1 exact match
```

### C. No-scale-up path

Choose Scale so `scale <= 1`.

Verify no unintended scalar multiplication or Scale change.

### D. Poisoned dormant backing rows

Start from a ciphertext that originally had high-Level backing storage, poison `r2...rL`, shrink it to Level 0 while preserving capacity, then run Fast ModUp basis.

Verify:

```text
q0/q1 match reference
poisoned high rows do not influence q0/q1
high rows are not deliberately synchronized
```

### E. Backing-storage reuse

Verify the common Bootstrap path can restore MaxLevel structural rows without allocating new O(N*L) coefficient backing when the Level-0 input came from a previously high-Level ciphertext.

A test may record underlying row pointers before shrink and confirm they are reused after ModUp restoration.

Also test fallback behavior for a genuinely Level-0 ciphertext that has insufficient capacity: correctness must still hold even if missing rows must be allocated.

### F. Validation

Reject at least:

```text
nil input
wrong degree
wrong ring type
input Level != 0
non-NTT input
Montgomery input
```

## Benchmarks

Benchmark Fast ModUp basis vs the Standard-style full-Q pre-Trace reference at:

```text
LogN = 13
LogN = 16
```

Report:

```text
ns/op
B/op
allocs/op
```

Benchmark the steady-state common path where high-Level backing already exists and is reused.

Setup/reset must be outside the timed region where appropriate and must not favor one implementation.

Also report one cold/fallback allocation observation for an input created genuinely at Level 0, but the primary Stage-A metric is steady-state Bootstrap reuse.

## Performance objective

Fast should avoid numerical work proportional to the full Q chain.

The expected hot numerical work is approximately:

```text
2 ciphertext components
x q0 INTT
x q0/q1 lift
x q0/q1 NTT
x optional q0/q1 scalar multiply
```

Structural high-Level rows may exist, but Fast must not NTT/INTT/arithmetic over q2...qL.

Do not require a fixed speedup threshold before measuring the implementation.

## Preserve existing behavior

Run:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/bootstrapping
go test ./circuits/ckks/dft
```

Do not regress:

```text
Fast ScaleDown
Fast Rescale
Fast evaluator
Fast DFT
Fast automorphism
```

## Out of scope

Do not implement in 006B:

```text
Trace
final public full ModUp semantics
DenseToSparse / SparseToDense
QP arithmetic
GadgetProduct
key switching
CoeffsToSlots orchestration
EvalMod
full Bootstrap
new sparse ciphertext representation
changes to ring.Poly Level semantics
```

## Documentation

Update the relevant status in `docs/FAST_CKKS_SPEC.md` to distinguish:

```text
Fast ModUp basis raise / scale alignment: implemented
Fast Trace: not implemented
full Fast ModUp: not implemented until Trace is integrated
```

Do not claim full Bootstrap support.

## Preferred files

Prefer a narrow change set around:

```text
circuits/ckks/bootstrapping/fast_modup.go
circuits/ckks/bootstrapping/fast_modup_test.go
circuits/ckks/bootstrapping/fast_scaledown.go   # only if FastEvaluator state needs a small extension
schemes/ckks/fast/evaluator.go                  # only if a narrow maintained-residue helper is needed
docs/FAST_CKKS_SPEC.md
```

Avoid unrelated core changes.

## Completion

Commit with:

```text
feat(ckks/fast): implement Fast ModUp basis raise
```

Push to:

```text
origin/fast-ckks
```

Final report should include:

1. Fast ModUp basis API/location;
2. q0-centered q1-lift strategy;
3. confirmation q2...qL are not numerically maintained;
4. backing-storage reuse strategy;
5. fallback behavior for fresh Level-0 storage;
6. scale-alignment behavior;
7. Standard/reference q0/q1 comparisons;
8. poisoned-residue result;
9. LogN13 benchmark;
10. LogN16 benchmark;
11. steady-state B/op and allocs/op;
12. cold/fallback allocation observation;
13. tests;
14. files changed;
15. commit hash;
16. push result.
