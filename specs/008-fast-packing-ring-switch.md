# Task 008 — Fast Packing / Unpacking + Ring-Degree Boundary

## Goal

Complete the Standard-ring Bootstrap boundary around `BootstrapMany` so the next task can wire the already-implemented Fast Bootstrap core end to end.

Current completed core:

```text
Fast ScaleDown          ✅
Fast ModUp + Trace      ✅
Fast CoeffsToSlots      ✅
Fast EvalMod            ✅
Fast SlotsToCoeffs      ✅
Fast ring-degree primitive ✅ (`FastN1ToN2` / `FastN2ToN1`)
```

Task 008 adds only the missing outer boundary:

```text
input ciphertexts in residual ring N1
    ↓
Fast Pack in N1 (when useful)
    ↓
Fast N1 → N2 ring-degree mapping (when N1 < N2)
    ↓
Fast Pack in N2
    ↓
[Bootstrap core — NOT wired here]
    ↓
Fast Unpack in N2
    ↓
Fast N2 → N1 mapping (when N1 < N2)
    ↓
Fast Unpack in N1
    ↓
original ciphertext count / LogSlots restored
```

Do **not** implement full `Bootstrap` / `BootstrapMany` orchestration in this task.

Base reviewed commit:

```text
302a18e89085ec3e256769df0c842a3c988e6aca
```

---

## 1. Source facts to preserve

Read current source before editing:

```text
circuits/ckks/bootstrapping/evaluator.go
schemes/ckks/fast/ring_degree.go
schemes/ckks/fast/ring_degree_test.go
circuits/ckks/bootstrapping/fast_scaledown.go
docs/FAST_CKKS_SPEC.md
```

Standard `BootstrapMany` currently does:

```text
PackAndSwitchN1ToN2
→ Evaluate each packed ciphertext
→ UnpackAndSwitchN2ToN1
→ restore residual DefaultScale
```

Standard packing itself performs staged sparse-polynomial combination:

```text
even += odd * X^k
```

using `MulCoeffsMontgomeryThenAdd`.

Standard unpacking creates shifted copies with inverse monomials `X^-k` and later restores the original `LogDimensions.Cols`.

Standard N1↔N2 switching currently calls `ApplyEvaluationKey`, which is not allowed in Fast.

The existing Fast primitive already provides the Stage-A zero-secret ring mapping:

```go
fastckks.FastN1ToN2(...)
fastckks.FastN2ToN1(...)
```

and explicitly reads/writes only q0/q1.

Do not reimplement this mathematics unless integration exposes a concrete defect.

---

## 2. Fast architecture contract

At this boundary:

```text
logical Level may be > 1
q0/q1 are authoritative
q2...qL may exist structurally but are dormant/stale
```

Therefore Fast packing/unpacking must never call full-level ring operations that implicitly consume every row.

For each monomial multiply/add, execute only on:

```text
component c0/c1
limb q0
limb q1
```

High rows must not influence q0/q1 and need not be synchronized.

No:

```text
Standard evaluator fallback
ApplyEvaluationKey
QP
GadgetProduct
key switching
CRT redistribution
per-coefficient big.Int
```

---

## 3. Required scope

Required Stage-A support:

```text
ring.Standard
Degree = 1
NTT-domain ciphertexts
q0/q1 storage present
consistent representation across the input slice
N2 == N1, or N2 == 2*N1
```

The `N2 == 2*N1` bound matches the already-tested Fast ring-degree primitive.

If `BootstrappingParameters.N() > 2*ResidualParameters.N()`, reject clearly in this task rather than silently chaining or falling back. Do not generalize ring-degree conversion unless an existing repository test/current target configuration proves it is necessary.

ConjugateInvariant is out of scope.

---

## 4. Fast packing monomial tables

Standard evaluator stores:

```text
xPow2N1
xPow2N2
xPow2InvN1
xPow2InvN2
```

Fast needs the same logical monomials but only for maintained q0/q1 arithmetic.

Prefer evaluator-owned precomputation/lazy initialization using q0/q1 only, conceptually from:

```go
params.RingQ().AtLevel(1)
```

while preserving the exact Standard exponent/index schedule.

Important: mirror Standard constructor arguments for N1/N2 forward and inverse tables. Do not derive a new packing convention.

Avoid regenerating monomial tables per call.

If minimal `NewFastEvaluator` constructions used by older ScaleDown/EvalMod tests do not provide `ResidualParameters`, packing-table initialization must not break them. Lazy packing initialization is acceptable/preferred.

---

## 5. Fast pack helper

Implement a Fast equivalent of Standard `pack`, preferably package-local under:

```text
circuits/ckks/bootstrapping/fast_packing.go
```

Conceptually:

```go
func (eval *FastEvaluator) fastPack(
    cts []rlwe.Ciphertext,
    ctxt packingContext,
    xPow2 []ring.Poly,
) ([]rlwe.Ciphertext, error)
```

Equivalent structure is acceptable.

Preserve Standard algorithm and metadata:

```text
logPackCTs = logMaxSlots - logSlots
logGap     = params.LogMaxSlots() - logSlots - 1

for each packing stage:
    even += odd * xPow2[logGap-i]

final packed ct.LogDimensions = ctxt.LogMaxDimensions
```

But the multiply-add must operate q0/q1 only.

Do not call:

```go
ctxt.Params.RingQ().AtLevel(level).MulCoeffsMontgomeryThenAdd(...)
```

on whole `ring.Poly` rows in production Fast packing.

Use per-maintained-limb subring operations or an equally explicit q0/q1 helper.

The representation (ordinary vs Montgomery ciphertext coefficients) must be preserved exactly; the monomial operand should be in the same Standard-compatible NTT/Montgomery convention as the current packing formula.

---

## 6. Fast unpack helper

Implement a Fast equivalent of Standard `unpack`.

Critical difference from Standard:

```go
ct.CopyNew()
```

must not be used to create every unpacked copy because it copies dormant high limbs.

Create independent structural ciphertexts, then copy only:

```text
q0/q1
metadata
```

For each required inverse-monomial multiplication, operate q0/q1 only.

Preserve Standard count/index logic:

```text
logPackCTs
n = min(NbPackedCTs, 1<<logPackCTs)
logGap
step/j/k traversal
```

Output ownership must be independent: changing/reusing one returned ciphertext must not corrupt its siblings.

---

## 7. Bootstrap-level packing APIs

Add Fast receiver methods mirroring Standard names/signatures where practical:

```go
func (eval *FastEvaluator) PackAndSwitchN1ToN2(
    cts []rlwe.Ciphertext,
) ([]rlwe.Ciphertext, *packingContext, *packingContext, error)

func (eval *FastEvaluator) UnpackAndSwitchN2ToN1(
    cts []rlwe.Ciphertext,
    ctxtN1, ctxtN2 *packingContext,
) ([]rlwe.Ciphertext, error)
```

### N1 == N2

No ring-degree conversion. Pack/unpack only in N2/current ring as Standard does.

### N2 == 2*N1

Preserve Standard ordering:

```text
Pack N1
→ FastN1ToN2
→ Pack N2
```

and reverse:

```text
Unpack N2
→ FastN2ToN1
→ Unpack N1
```

Allocate structural destination ciphertexts for ring-degree conversion, but do not copy dormant q2+ values into them.

Use `ResidualParameters.RingQ().AtLevel(level)` and `BootstrappingParameters.RingQ().AtLevel(level)` consistently with the current ciphertext level.

Do not invoke evaluation keys: the current Fast zero-secret Stage-A semantics use the pure ring-degree mapping already implemented.

---

## 8. Metadata contract

Preserve at least:

```text
Scale
IsNTT
IsMontgomery
IsBatched
IsBitReversed
LogDimensions
Degree
Level
```

Packing must set `LogDimensions` exactly as Standard packing does.

Final unpack must restore:

```go
ctsOut[i].LogDimensions.Cols = original LogSlots
```

Do not invent scale changes in packing/ring-switching. These are representation/layout operations, not CKKS rescaling operations.

---

## 9. Validation

Reject clearly:

```text
nil/empty ciphertext slice
nil metadata
ring type != Standard
Degree != 1
non-NTT input
missing q0/q1 storage
inconsistent N across slice
inconsistent LogSlots across slice
inconsistent Level / domain representation when required
LogSlots > target LogMaxDimensions.Cols
N2 < N1
N2 > 2*N1 for current Stage-A scope
q0/q1 moduli mismatch across N1/N2 at the used Level
invalid/missing packing context during unpack
```

Do not mutate the branch into a generic packing library.

---

## 10. Correctness tests

### A. Fast pack exact q0/q1 oracle

For valid sparse ciphertext inputs in one ring, compare Fast packing against the current Standard `pack` mathematics or an exact ring reference.

Compare exact:

```text
q0/q1
Level
Scale
LogDimensions
representation flags
```

Use at least:

```text
2 ciphertexts
3 ciphertexts (odd count)
```

### B. Fast unpack oracle

Compare Fast unpack against Standard/reference for q0/q1 and metadata.

Exercise more than one unpack stage.

### C. poison test

Use two inputs with identical q0/q1 but different q2...qL poison.

Verify complete Fast pack output q0/q1 is identical.

Repeat for unpack.

### D. same-ring pack→unpack

Use valid sparse inputs and verify the packing boundary restores the same logical set/reference values and original `LogSlots`.

Do not use arbitrary dense random ciphertexts if the Standard packing algorithm assumes sparse-layout inputs.

### E. N1→N2 boundary

With `N2 = 2*N1`:

```text
Fast Pack N1
→ FastN1ToN2
→ Fast Pack N2
```

verify q0/q1 against a source-faithful reference assembled from Standard/raw ring operations.

### F. reverse N2→N1 boundary

Verify unpack/switch/unpack restores expected q0/q1 and metadata.

### G. counts

At minimum cover:

```text
1
2
3
5
```

ciphertexts where parameter dimensions permit.

### H. returned ownership

Run packing/unpacking repeatedly and verify earlier returned outputs are not mutated by later calls.

---

## 11. Existing ring-degree primitive

Re-run existing tests for:

```text
FastN1ToN2
FastN2ToN1
NTT path
poisoned dormant limbs
```

Do not rewrite it solely for style.

If integration exposes repeated allocation in `ring_degree.go` from per-call temporary q0/q1 polys, report benchmark impact first. A small evaluator-owned scratch integration is allowed only if clearly local and correctness-preserving.

Do not open an optimization subtask automatically.

---

## 12. Benchmarks

Benchmark warm packing boundary for representative counts, preferably:

```text
LogN13, same ring, 2 or 4 input ciphertexts
LogN13, N1→N2 factor-2 case if practical
```

If LogN16 is practical, include it; otherwise prioritize functional integration because 009 will benchmark the end-to-end path.

Report:

```text
ns/op
B/op
allocs/op
number of input ciphertexts
N1/N2
```

Do not claim a speedup over Standard unless the compared benchmark measures the same mathematical/layout operation.

---

## 13. Regression commands

Run at least:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/polynomial
go test ./circuits/ckks/mod1
go test ./circuits/ckks/dft
go test ./circuits/ckks/bootstrapping
```

No regressions to 006/007 stages.

---

## 14. Out of scope

Do not implement:

```text
Fast Bootstrap / BootstrapMany orchestration
ConjugateInvariant packing
RealToComplex / ComplexToReal
DenseToSparse / SparseToDense
arbitrary N2/N1 > 2 ratio
Standard evaluation-key switching
QP/GadgetProduct
new ciphertext storage representation
noise-fidelity work
final performance cleanup
```

Task 009 will connect:

```text
Pack/switch
→ Fast bootstrap core
→ unpack/switch
```

---

## Preferred files

Prefer:

```text
circuits/ckks/bootstrapping/fast_packing.go
circuits/ckks/bootstrapping/fast_packing_test.go
circuits/ckks/bootstrapping/fast_packing_bench_test.go
docs/FAST_CKKS_SPEC.md
```

Small changes to `fast_scaledown.go` / Fast evaluator fields/constructor are acceptable for reusable packing tables.

Avoid changing Standard `evaluator.go`.

Avoid changing `schemes/ckks/fast/ring_degree.go` unless integration proves a concrete need.

---

## Completion

Commit with:

```text
feat(ckks/fast): add Fast bootstrap packing boundary
```

Push to:

```text
origin/fast-ckks
```

Final report must include:

1. Fast packing API/files;
2. q0/q1 monomial-table strategy;
3. Fast pack maintained-limb operation;
4. Fast unpack copy/allocation strategy;
5. N1==N2 behavior;
6. N2==2*N1 ring-switch integration;
7. confirmation no evaluation-key/QP/Gadget fallback;
8. metadata behavior;
9. exact pack oracle results;
10. exact unpack oracle results;
11. poison test results;
12. same-ring round-trip result;
13. factor-2 N1/N2 boundary result;
14. odd-count/multi-count results;
15. ownership/reuse result;
16. benchmark results;
17. remaining structural allocation limitations;
18. regression tests;
19. files changed;
20. commit hash;
21. push result.
