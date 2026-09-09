# Task 009-P — End-to-End Fast Bootstrap Review Fix

## Goal

Close the remaining verification gaps in Task 009 without redesigning the accepted production Fast Bootstrap path.

Reviewed implementation commit:

```text
79649cfb2af751b41f92358f8330cbbf2c6996a3
feat(ckks/fast): wire end-to-end Fast Bootstrap
```

Independent source review did not identify a production mathematical defect in the current orchestration. The blocking issue is incomplete acceptance evidence for several end-to-end branches required by the 009 spec.

Unless a stronger test exposes a real bug, prefer test/benchmark-only changes.

---

## Finding 1 — no full-slot real+imag end-to-end Bootstrap test

Current `fast_bootstrap_test.go` uses `LogN=4` with `LogSlots=2` or `1`.

For a Standard ring with `LogN=4`, maximum complex slots are `LogSlots=3`, so these tests exercise the sparse/repacked path. They do not exercise the branch where Fast `CoeffsToSlotsNew` returns both:

```text
ctReal != nil
ctImag != nil
```

### Required fix

Add a full-slot end-to-end case with:

```text
LogSlots == BootstrappingParameters.LogMaxSlots()
```

Use deterministic nonzero complex values with nonzero real and imaginary parts.

Verify:

```text
C2S produces ctReal and ctImag
EvalMod is applied to both
S2C recombines them
Bootstrap succeeds
final output IsNTT = true
final output IsMontgomery = false
final Level = ResidualParameters.MaxLevel()
final Scale = ResidualParameters.DefaultScale()
decoded output ~= original values
```

Tolerance should be source/test justified; `1e-2` is acceptable if the observed errors support it.

Do not mock the branch. Exercise the actual public `Bootstrap` path.

---

## Finding 2 — Level-1 input test is not numerical

Current:

```text
TestFastBootstrapLevelOneInput
```

checks output Level and input ownership, but does not decode the result.

### Required fix

Use deterministic nonzero complex values at Level 1 and verify:

```text
Fast Bootstrap decoded output ~= original values
```

with documented tolerance.

Keep the existing ownership assertion.

This test must actually exercise the Level-1 → ScaleDown → Level-0 transition rather than converting the source to Level 0 before Bootstrap.

---

## Finding 3 — BootstrapMany odd-count test does not verify message order

Current 3-ciphertext test checks:

```text
count
LogSlots
Level
IsMontgomery
```

but does not verify that each output still corresponds to the correct input after packing/unpacking.

### Required fix

Use 3 distinct deterministic message vectors.

After `BootstrapMany`, decode every output and compare it to its corresponding original vector.

This must catch:

```text
reordering
cross-ciphertext mixing
wrong odd-tail unpacking
aliasing between outputs
```

Also verify returned ciphertext ownership by mutating one output after the call and confirming the others remain unchanged.

---

## Finding 4 — prove compatibility with ordinary decrypt/decode API

Current Fast numerical helper directly treats `c0` as plaintext, which is valid as a zero-secret semantic oracle, but the 009 public-output contract also requires evidence that the returned ordinary/non-Montgomery ciphertext can pass through the normal decrypt/decode API.

### Required fix

Construct an explicit zero secret key under the residual parameters and use:

```go
rlwe.NewDecryptor(residualParams, zeroSK).DecryptNew(fastOut)
ckks.NewEncoder(residualParams).Decode(...)
```

Compare to the original deterministic values.

Do not use a nonzero secret for Fast Stage-A output unless the current production semantics support it.

This test is evidence for the public representation boundary:

```text
NTT = true
Montgomery = false
Level 0 or 1
maintained public residues ordinary
```

---

## Finding 5 — representative benchmark is not the requested default-like Mod1 profile

Current benchmark helper uses approximately:

```text
Mod1Degree = 30
K = 4
DoubleAngle = 2
```

The 009 spec requested, when practical, a default-like profile around:

```text
Mod1Degree = 30
K = 16
DoubleAngle = 3
Mod1InvDegree = 0
```

### Required fix

Add or change the end-to-end benchmark to use:

```text
CosDiscrete
K = 16
Mod1Degree = 30
DoubleAngle = 3
Mod1InvDegree = 0
```

with enough Q levels for the actual DFT + EvalMod + S2C depth.

Run warm:

```text
LogN13
LogN16
```

Report:

```text
ns/op
B/op
allocs/op
LogN
LogSlots
input Level
output Level
```

Use normal Go benchmark timing if practical. If `-benchtime=1x` is necessary due runtime, state that clearly and do not overinterpret small timing deltas.

This is still 009 correctness/integration benchmarking. Task 010 will do final profiling and optimization.

---

## DFT planning parity evidence

Task 009 correctly extracted the Standard planning/matrix generation into shared:

```text
buildBootstrapCircuitData
```

and Standard `initialize` now calls the same helper.

Add one focused test confirming a Standard evaluator and Fast evaluator initialized from the same `Parameters` expose equal source-level planning data for at least:

```text
Mod1Parameters.LevelQ
Mod1Parameters.QDiff
C2S matrix LevelQ / level schedule
S2C matrix LevelQ / level schedule
adjusted C2S scaling
adjusted S2C scaling
```

Do not duplicate the formulas in the test if direct state comparison is possible.

---

## Factor-two N1→N2 path

Production code claims support for:

```text
N2 = 2*N1
```

Task 008 already provides exact q0/q1 forward/reverse packing-boundary references, and Task 009 added Level-0 ring-degree tests.

Preferred if practical: add one small end-to-end factor-two `Bootstrap` or `BootstrapMany` numerical smoke test.

If runtime/setup makes that disproportionate, it is acceptable to retain the source-backed combination of:

```text
Task-008 exact factor-two pack/switch/unpack tests
+
Task-009 Level-0 ring-degree tests
+
same-ring full Bootstrap numerical tests
```

but the final report must explicitly state that factor-two was not numerically exercised through the entire Bootstrap core.

Do not silently claim a full factor-two numerical oracle if none exists.

---

## Preserve accepted production architecture

Unless a new test fails, do not change:

```text
buildBootstrapCircuitData
bootstrapCore
Bootstrap
BootstrapMany
public IMForm finalization
Fast DFT execution
Fast EvalMod execution
Fast packing/ring-degree execution
```

No Standard evaluator, QP, GadgetProduct, evaluation-key, or key-switch fallback may be added to production Fast execution.

---

## Required regressions

Run at least:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/polynomial
go test ./circuits/ckks/mod1
go test ./circuits/ckks/dft
go test ./circuits/ckks/bootstrapping
```

Run the new full-slot, Level-1 numerical, BootstrapMany-order, zero-secret decryptor, and benchmark cases.

---

## Preferred files

Prefer changes limited to:

```text
circuits/ckks/bootstrapping/fast_bootstrap_test.go
circuits/ckks/bootstrapping/fast_bootstrap_bench_test.go
```

A tiny production fix is allowed only if a stronger oracle exposes a real defect.

---

## Completion

Commit with:

```text
test(ckks/fast): close end-to-end Bootstrap review gaps
```

Push to:

```text
origin/fast-ckks
```

Final report must include:

1. full-slot real+imag end-to-end numerical result;
2. Level-1 numerical result;
3. BootstrapMany 3-input order/message result;
4. output ownership result;
5. ordinary zero-secret decryptor + encoder decode result;
6. shared DFT planning parity result;
7. factor-two full-core numerical result, or explicit statement that coverage remains compositional;
8. default-like K16/degree30/DoubleAngle3 benchmark parameters;
9. LogN13 benchmark;
10. LogN16 benchmark;
11. B/op and allocs/op;
12. whether production code changed;
13. regressions;
14. files changed;
15. commit hash;
16. push result.
