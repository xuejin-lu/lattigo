# Task 007B-P — Fast EvalMod Review Fix

## Goal

Close the remaining verification gaps in Task 007B without redesigning the accepted production Fast Mod1 path.

Reviewed implementation commit:

```text
de606e11316bcb1f29f1fa606adab078021d8ae8
feat(ckks/fast): implement Fast EvalMod
```

The production implementation currently mirrors the Standard Mod1 flow closely enough that no source-backed mathematical defect has been identified. The blocking issue is insufficient correctness evidence for the full supported path.

Do not begin full Bootstrap orchestration.

---

## Blocking finding 1 — `MatchesStandard` does not numerically compare Standard and Fast

Current:

```text
TestFastMod1MatchesStandardCosDiscrete
```

constructs a Standard Mod1 result and a Fast Mod1 result, but only checks:

```text
Level
Scale
Fast NTT/Montgomery flags
```

It does **not** decode and compare Standard vs Fast values.

The only numerical assertion is:

```text
zero message -> Fast output approximately zero
```

This is too weak for EvalMod because an implementation can have a wrong normalization, offset, coefficient scale, or DoubleAngle behavior and still map the all-zero test case near zero.

### Required fix

Use nonzero deterministic CKKS slot values in the supported input domain.

Run the same valid full-Q input through:

```text
Standard Mod1 EvaluateNew
Fast Mod1 EvaluateNew
```

Decode both outputs appropriately and compare:

```text
Fast decoded slots ~= Standard decoded slots
```

with a documented CKKS tolerance.

Do not compare only metadata.

The input must be nontrivial enough that omission of normalization, cosine offset, polynomial evaluation, or DoubleAngle would change the expected values.

---

## Blocking finding 2 — no default-like degree-30 / DoubleAngle-3 numerical test

The current default-like configuration:

```text
K = 16
Mod1Degree = 30
DoubleAngle = 3
Mod1InvDegree = 0
```

is used in benchmark code only.

Functional numerical testing currently uses a much smaller configuration such as degree 6 / DoubleAngle 2.

### Required fix

Add one focused default-like numerical correctness test using:

```text
CosDiscrete
K = 16
Mod1Degree = 30
DoubleAngle = 3
Mod1InvDegree = 0
```

with sufficient Q levels.

Preferred oracle:

```text
Standard Mod1 decoded result
vs
Fast Mod1 decoded result
```

on the same deterministic nonzero zero-secret-style input.

Verify at least:

```text
output Level matches Standard
output Scale matches Standard before Bootstrap wrapper
decoded values within tolerance
Degree = 1
NTT/Montgomery preserved
```

Then test the Bootstrap Fast `EvalMod` wrapper separately for:

```text
Scale = BootstrappingParameters.DefaultScale()
```

Do not use the benchmark itself as correctness evidence.

---

## Blocking finding 3 — direct oracle for new Fast `MulThenAdd` scale promotion

Task 007B modified:

```text
schemes/ckks/fast.Evaluator.MulThenAdd
```

so that when the accumulator scale is lower than the term scale, the accumulator is promoted with q0/q1-only integer scaling before the fused add.

This behavior is source-motivated and consistent with the Standard scale-alignment pattern, but the new branch needs a direct regression test.

### Required test

Construct a scalar `MulThenAdd` case where:

```text
op0.Scale > opOut.Scale
```

and where the scale-promotion branch definitely executes.

Use the same q0/q1 inputs with a Standard CKKS evaluator and Fast evaluator.

Compare:

```text
exact q0/q1 residues
output Scale
output Level
```

after accounting for representation consistently.

Also poison q2...qL in Fast inputs/output and verify the q0/q1 answer is unchanged.

Test both if practical:

```text
integer scalar
non-integer scalar
```

At minimum the non-integer scalar case used by polynomial coefficient accumulation must be covered.

Do not broaden this into a generic evaluator redesign.

---

## Normalization / cosine-offset sensitivity

The original 007B spec required evidence that the normalization and cosine offset are not dead code.

Add either:

### Preferred
A numerical test with nonzero inputs where Standard vs Fast comparison would fail if the Fast normalization or cosine offset were intentionally omitted.

### Or
A small dedicated source-faithful reference test that independently computes the prepared input state after:

```text
res.Scale = ScalingFactor()
Fast Add(cosine offset)
```

and compares exact q0/q1 maintained residues + metadata before polynomial evaluation.

Do not duplicate the whole Mod1 implementation as the reference.

---

## Preserve accepted production architecture

Unless the stronger tests reveal a real discrepancy, do not change:

```text
circuits/ckks/mod1/fast.go
```

except for small testability fixes.

Keep:

```text
Standard ring
CosDiscrete only
Mod1InvPoly == nil
scaling = 1
exact LevelQ input
NTT + Montgomery
q0/q1 authoritative
q2...qL dormant
Fast polynomial evaluator reuse
Fast DoubleAngle primitives
final Mod1 input-scale restoration
Bootstrap wrapper default-scale restoration
```

No Standard/QP/Gadget/key-switch fallback may be added to production.

---

## Required regression tests

Run at least:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/polynomial
go test ./circuits/ckks/mod1
go test ./circuits/ckks/dft
go test ./circuits/ckks/bootstrapping
```

Re-run the existing Fast EvalMod benchmarks:

```text
LogN13
LogN16
```

Report `ns/op`, `B/op`, and `allocs/op`.

No performance improvement is required in 007B-P unless the correctness fix itself changes the hot path.

---

## Preferred files

Prefer changes limited to:

```text
circuits/ckks/mod1/fast_test.go
circuits/ckks/bootstrapping/fast_evalmod_test.go
schemes/ckks/fast/evaluator_ntt_test.go
```

Modify production files only if the stronger oracle exposes a real bug.

Do not modify Standard behavior.

---

## Completion

Commit with:

```text
test(ckks/fast): strengthen Fast EvalMod oracles
```

Push to:

```text
origin/fast-ckks
```

Final report must include:

1. nonzero Standard-vs-Fast Mod1 numerical comparison;
2. test input construction and tolerance;
3. default-like degree-30 / DoubleAngle-3 numerical result;
4. output Level/Scale comparison;
5. normalization/cosine-offset sensitivity evidence;
6. direct `MulThenAdd` scale-promotion oracle;
7. poisoned q2...qL result for scale promotion;
8. whether production code required changes;
9. LogN13 benchmark;
10. LogN16 benchmark;
11. B/op and allocs/op;
12. regression tests;
13. files changed;
14. commit hash;
15. push result.
