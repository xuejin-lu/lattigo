# Task 006C — Implement Fast Trace and Complete Fast ModUp Boundary

## Goal

Complete the Fast Bootstrap `ModUp` boundary by adding Fast Trace and the final representation handoff required by Fast DFT.

Starting point:

```text
Fast ScaleDown
-> Level 0
-> NTT
-> ordinary/non-Montgomery
```

Task 006B already provides the private pre-Trace basis raise:

```text
modUpBasis
-> MaxLevel structure
-> q0/q1 numerically maintained
-> q2...qL dormant/stale
-> NTT
-> ordinary/non-Montgomery
```

Task 006C must produce the complete Fast ModUp boundary:

```text
ScaleDown output
-> modUpBasis
-> Fast Trace
-> q0/q1 Montgomery conversion
-> MaxLevel, NTT, Montgomery
-> ready for Fast CoeffsToSlots
```

Base commit:

```text
9b69fef1ffdc20a9befc050bbe6e61c6a4ebde3b
feat(ckks/fast): implement Fast ModUp basis raise
```

Do not implement EvalMod or full Bootstrap orchestration in this task.

---

## 1. Implement Fast Trace in `schemes/ckks/fast`

Add a bounded Fast Trace operation, preferably:

```go
func (eval *Evaluator) Trace(ctIn *rlwe.Ciphertext, logN int, opOut *rlwe.Ciphertext) error
```

It must implement the same mathematical Trace as `core/rlwe/inner_sum.go`, but only on maintained q0/q1 residues and without evaluation keys, GadgetProduct, QP arithmetic, or Standard evaluator fallback.

Primary supported contract:

```text
ring.Standard
Degree = 1
Level >= 1
NTT domain
q0/q1 valid
q2...qL dormant
```

Trace may support either ordinary or Montgomery NTT representation as long as input/output representation flags match and q0/q1 semantics remain correct.

Do not add non-NTT support unless it is essentially free and does not introduce full-Q work. The Bootstrap path entering Trace is NTT.

---

## 2. Preserve Standard Trace semantics exactly

Read the current source in:

```text
core/rlwe/inner_sum.go
```

The Standard logic is:

```go
gap := 1 << (params.LogN() - logN - 1)
if logN == 0 {
    gap <<= 1
}
```

If `gap <= 1`, Trace is a copy/no-op.

If `gap > 1`:

1. Pre-multiply by `gap^-1` modulo the current ciphertext modulus.
2. For:

   ```go
   i := logN; i < params.LogN()-1; i++
   ```

   apply:

   ```text
   ct <- ct + phi_{GaloisElement(1<<i)}(ct)
   ```

3. If:

   ```text
   logN == 0 && ring.Standard
   ```

   apply the final order-two automorphism and add once more.

Metadata Scale must remain unchanged by Trace.

Validate `logN` so invalid values return an error instead of relying on invalid shifts.

---

## 3. Fast inverse scaling

The Standard implementation computes an integer inverse of `gap` modulo `Q_level` and multiplies the ciphertext by it.

Fast may compute the same one-time `big.Int` inverse using:

```text
RingQ.ModulusAtLevel[level]
```

and pass it through the existing q0/q1-only integer multiplication path.

This is acceptable because it is one scalar operation per Trace, not per-coefficient arbitrary-precision CRT.

Do not use per-coefficient `big.Int`.

The resulting q0 and q1 residues must equal multiplication by `gap^-1 mod qi` for `i in {0,1}`.

Do not change ciphertext Scale metadata during this normalization.

---

## 4. Fuse automorphism + add without full ciphertext scratch

Do not allocate a MaxLevel ciphertext for every Trace step.

The Fast evaluator already owns:

```text
automorphism index cache
q0/q1 automorphism scratch
q0/q1 NTT scratch
```

Use evaluator-owned q0/q1 scratch to implement each Trace step conceptually as:

```text
for each ciphertext component:
    tmp = phi(component)
    component = component + tmp  mod q0/q1
```

Only q0/q1 may be read or written.

`q2...qL` must remain untouched/dormant.

A first call may populate the existing automorphism index cache. Repeated Trace calls with the same Galois elements must reuse the cache instead of rebuilding O(N) index arrays every Bootstrap.

Do not add a global mutable cache.

---

## 5. Aliasing and output storage

Support:

```text
in-place Trace
out-of-place Trace
```

For out-of-place operation, require structurally compatible output storage and copy only the maintained q0/q1 state plus metadata. Do not synchronize dormant high residues merely for raw ciphertext equality.

Preserve:

```text
Level
Scale
IsNTT
IsMontgomery
LogDimensions / metadata
```

---

## 6. Complete public Fast `ModUp`

After Fast Trace is implemented, add the public Bootstrap boundary:

```go
func (eval *FastEvaluator) ModUp(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error)
```

Its sequence must be:

```text
modUpBasis(ct)
-> FastCKKS.Trace(ct, CoeffsToSlotsParameters.LogSlots, ct)
-> convert maintained q0/q1 to Montgomery form
-> return ct
```

Do not call Standard Bootstrap `ModUp`.

Do not perform DenseToSparse / SparseToDense switching.

The current Stage-A Fast configuration remains the no-key ModUp branch.

---

## 7. Fast DFT handoff requirement

Current Fast DFT explicitly requires:

```text
NTT
Montgomery
Degree 1
Level >= 1
```

Therefore the complete Fast `ModUp` output must be:

```text
Level = BootstrappingParameters.MaxLevel()
IsNTT = true
IsMontgomery = true
q0/q1 valid in Montgomery representation
q2...qL dormant/stale
```

The Montgomery conversion must touch only q0/q1.

Do not MForm q2...qL.

A narrow private helper in Bootstrap or a narrowly-scoped maintained-residue helper in `schemes/ckks/fast` is acceptable. Do not broaden unrelated evaluator APIs.

`modUpBasis` itself may remain ordinary/non-Montgomery so its 006B contract and tests stay explicit.

---

## 8. Trace reference oracle

Do not compare raw Fast ciphertexts against Standard `rlwe.Evaluator.Trace` ciphertext components, because Standard automorphism performs key switching and therefore raw c0/c1 are not an appropriate equality oracle.

Instead add a test-only ring-level mathematical Trace reference that:

```text
uses the same gap and Galois elements as Standard Trace
applies direct ring automorphisms
adds the results
may operate full-Q or q0/q1 in test/reference code
```

Compare exact maintained q0/q1 residues against this reference.

The reference must be derived from the current Standard source, not invented from memory.

---

## 9. Required Trace tests

Add focused tests for at least:

### A. `gap == 1`

Use max-slot-like `logN` where Trace is a no-op.

Verify in-place and out-of-place q0/q1 equality and unchanged Scale.

### B. Nontrivial Trace

Choose `logN` so multiple automorphism/add stages execute.

Verify exact q0/q1 equality against the ring-level reference.

### C. `logN == 0`

Exercise the special final order-two automorphism.

Verify exact q0/q1 equality against reference.

### D. Inverse normalization

Use nontrivial q0/q1 values and verify the `gap^-1` pre-multiplication is present. A test that would fail if normalization were omitted is required.

### E. Poisoned dormant residues

Poison q2...qL before Trace.

Verify:

```text
q0/q1 match reference
q2...qL do not influence q0/q1
high rows are not deliberately synchronized
```

### F. Representation

Exercise ordinary NTT Trace. If Montgomery Trace is supported, test it too.

### G. Validation

Reject at least:

```text
nil evaluator/input/output
wrong degree
wrong ring type
Level < 1
non-NTT input if non-NTT is intentionally unsupported
mismatched input/output Level or representation
invalid logN
```

---

## 10. Required full Fast ModUp tests

Test the new public `FastEvaluator.ModUp` from a valid Level-0 ordinary NTT input.

Verify:

```text
modUpBasis semantics preserved
Trace semantics preserved
Level = MaxLevel
Scale matches Standard-style ModUp basis + Trace semantics
IsNTT = true
IsMontgomery = true
q0/q1 equal the Montgomery conversion of the mathematical reference
q2...qL remain dormant
```

Test both:

```text
Trace no-op / gap==1
nontrivial Trace / gap>1
```

Also verify that backing-row reuse from 006B still works through the complete ModUp call.

Do not require q2...qL equality with Standard.

---

## 11. Performance tests

Benchmark at:

```text
LogN = 13
LogN = 16
```

### Trace benchmark

Use a nontrivial `logN` that executes multiple Trace automorphisms.

Report warm steady-state Fast Trace:

```text
ns/op
B/op
allocs/op
```

Warm the automorphism index cache before timing so the primary benchmark represents repeated Bootstrap execution.

Also report the first/cold call separately if useful; do not mix one-time index-cache creation into the steady-state number.

### Complete ModUp benchmark

Benchmark:

```text
modUpBasis + Trace + q0/q1 Montgomery conversion
```

for:

```text
LogN13
LogN16
```

Use the steady-state backing-storage reuse path from 006B.

Compare against an appropriate full-Q mathematical/reference implementation when practical, but do not construct Standard evaluation keys solely to manufacture a benchmark. Clearly label reference benchmarks as reference rather than actual Standard key-switched ModUp if that is what is measured.

The main Stage-A requirement is that Fast Trace itself performs no q2...qL arithmetic and no QP/key-switch work.

---

## 12. Preserve existing tests

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
Fast ModUp basis
Fast automorphism
Fast DFT
Fast LinearTransform
```

---

## 13. Out of scope

Do not implement in 006C:

```text
EvalMod
polynomial PowerBasis integration
CoeffsToSlots orchestration inside Bootstrap
SlotsToCoeffs orchestration
full Bootstrap Evaluate/Bootstrap
DenseToSparse / SparseToDense
QP arithmetic
GadgetProduct
Standard key switching
new sparse ciphertext representation
ring.Poly / rlwe.Element Level redesign
```

Do not modify Standard behavior.

---

## 14. Documentation

Update `docs/FAST_CKKS_SPEC.md` so status becomes conceptually:

```text
Fast ScaleDown: implemented
Fast ModUp basis raise: implemented
Fast Trace: implemented
Fast ModUp boundary: implemented for Stage-A no-key Standard-ring path
Fast CoeffsToSlots adapter: already implemented, not yet fully orchestrated
Fast EvalMod: not implemented
full Fast Bootstrap: not implemented
```

Do not claim full Bootstrap completion.

---

## Preferred files

Prefer a narrow change set around:

```text
schemes/ckks/fast/trace.go
schemes/ckks/fast/trace_test.go
circuits/ckks/bootstrapping/fast_modup.go
circuits/ckks/bootstrapping/fast_modup_test.go
docs/FAST_CKKS_SPEC.md
```

Small changes to existing Fast evaluator scratch/helpers are allowed if required for an allocation-free Trace hot path.

Avoid unrelated cleanup.

---

## Completion

Commit with:

```text
feat(ckks/fast): implement Fast Trace and ModUp
```

Push to:

```text
origin/fast-ckks
```

Final report should include:

1. Fast Trace API/location;
2. exact `gap` / inverse normalization behavior;
3. Galois elements used;
4. handling of `logN == 0`;
5. automorphism scratch/cache reuse strategy;
6. confirmation q2...qL are untouched;
7. in-place/out-of-place result;
8. public Fast ModUp sequence;
9. Montgomery handoff strategy;
10. q0/q1 reference comparisons;
11. poisoned-residue result;
12. backing-row reuse result;
13. Trace LogN13 benchmark;
14. Trace LogN16 benchmark;
15. complete ModUp LogN13 benchmark;
16. complete ModUp LogN16 benchmark;
17. B/op and allocs/op;
18. tests;
19. files changed;
20. commit hash;
21. push result.
