# Task 007B — Fast EvalMod / Mod1

## Goal

Implement the Stage-A Fast homomorphic modular-reduction boundary used by CKKS Bootstrap:

```text
Fast CoeffsToSlots output
    -> Fast Mod1 / EvalMod
    -> ready for Fast SlotsToCoeffs
```

Base implementation commit:

```text
b0a43f98a019ecd61e13a987727466f68fdf6e8e
test(ckks/fast): close Fast polynomial review gaps
```

The Fast polynomial execution engine from 007A/007A-P is accepted and must be reused.

Task 007B must preserve the current Standard Mod1 mathematics while replacing Standard ciphertext execution with the explicit q0/q1 Fast arithmetic path.

Do not implement full Bootstrap orchestration in this task.

---

## 1. Source of truth

Before editing, read current source at least in:

```text
circuits/ckks/mod1/mod1_evaluator.go
circuits/ckks/mod1/mod1_parameters.go
circuits/ckks/polynomial/fast.go
circuits/ckks/bootstrapping/evaluator.go
circuits/ckks/bootstrapping/fast_scaledown.go
circuits/ckks/bootstrapping/fast_modup.go
docs/FAST_CKKS_SPEC.md
```

Standard Mod1 `EvaluateAndScaleNew` currently performs, conceptually:

```text
1. normalize ciphertext interpretation to mod 1 by changing Scale metadata
2. compute target polynomial scale for the later DoubleAngle chain
3. apply the cosine offset / Chebyshev change-of-variable preparation
4. evaluate Mod1Poly
5. run DoubleAngle repeatedly:
       y <- y^2
       y <- 2y
       y <- y - sqrt2pi
       Rescale
6. optional inverse/ArcSine polynomial
7. restore the input/output Scale interpretation
```

Standard Bootstrap `EvalMod` then resets the returned Scale to:

```text
BootstrappingParameters.DefaultScale()
```

Fast must preserve these semantics for the bounded Stage-A configuration.

---

## 2. Bounded Stage-A scope

Required in 007B:

```text
ring.Standard
Mod1Type = mod1.CosDiscrete
Mod1InvPoly == nil / Mod1InvDegree = 0
single EvaluateNew / scaling = 1 path
Degree-1 ciphertext
Level = Mod1Parameters.LevelQ at the Fast Mod1 entry
NTT
Montgomery
q0/q1 authoritative
q2...qL dormant
```

The project target uses the cosine path with polynomial evaluation followed by DoubleAngle; the current engineering profile uses degree ~30 and DoubleAngle = 3.

Out of scope for 007B:

```text
SinContinuous
CosContinuous
Mod1InvPoly / ArcSine
EvaluateAndScaleNew with scaling != 1
iterative/META bootstrapping scaling
PolynomialVector
full Bootstrap orchestration
```

Reject unsupported configurations explicitly. Do not silently fall back to Standard.

---

## 3. Instantiate the complete Mod1 parameters in Fast Bootstrap

Current `bootstrapping.NewFastEvaluator` only constructs a partial `mod1.Parameters` containing the fields needed by ScaleDown/ModUp.

For 007B, construct the complete Mod1 parameter object from the current source-of-truth literal:

```go
mod1.NewParametersFromLiteral(
    params.BootstrappingParameters,
    params.Mod1ParametersLiteral,
)
```

or the exact equivalent supported by current source.

This complete object must provide the actual:

```text
LevelQ
LogDefaultScale
Mod1Type
LogMessageRatio
DoubleAngle
QDiff
Sqrt2Pi
Mod1Poly
Mod1InvPoly
K / interval
```

ScaleDown and ModUp must continue using the same resulting `ScalingFactor()` and `MessageRatio()` semantics.

Do not duplicate Mod1 parameter-generation formulas in Bootstrap.

---

## 4. Add a bounded Fast Mod1 evaluator

Prefer adding:

```text
circuits/ckks/mod1/fast.go
```

with an API conceptually like:

```go
type FastEvaluator struct {
    Parameters          Parameters
    FastCKKS            *fastckks.Evaluator
    PolynomialEvaluator *ckkspolynomial.FastEvaluator
}

func NewFastEvaluator(
    eval *fastckks.Evaluator,
    evalPoly *ckkspolynomial.FastEvaluator,
    params Parameters,
) *FastEvaluator

func (eval *FastEvaluator) EvaluateNew(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error)
```

Equivalent narrow naming is acceptable.

This evaluator must not embed or instantiate a Standard `ckks.Evaluator`.

It must not implement `schemes.Evaluator`.

---

## 5. Preserve the exact Standard normalization / change-of-variable sequence

Do not reinterpret the Mod1 mathematics from memory. Mirror current `mod1.Evaluator.EvaluateAndScaleNew` source for `scaling = 1` and `Mod1InvPoly == nil`.

Important distinction:

The generic CKKS polynomial wrapper does **not** automatically apply Chebyshev interval change-of-basis. Its current source explicitly requires the caller to prepare the ciphertext before polynomial evaluation.

Therefore Fast Mod1 must preserve the same caller-side preparation as Standard, including:

### A. Scale normalization

Conceptually Standard changes:

```go
res.Scale = evm.ScalingFactor()
```

without changing polynomial coefficients.

This is a CKKS interpretation change, not a full-Q arithmetic operation.

Use the same metadata semantics.

### B. targetScale calculation

Compute the polynomial target scale using the same current-Q schedule and repeated square-root logic as Standard:

```text
for each DoubleAngle stage:
    targetScale *= corresponding Qi
    targetScale = sqrt(targetScale)
```

Use the exact current indices from Standard source; do not hardcode the default chain layout.

### C. cosine offset

For `CosDiscrete`, apply the exact current Standard offset derived from:

```text
Mod1Poly interval
IntervalShrinkFactor = 2^DoubleAngle
```

using Fast scalar Add.

Do not add another change-of-basis scalar unless current Standard source actually does so. The preparation is distributed across Bootstrap/Mod1 scale normalization and the offset; preserve source behavior rather than re-deriving an alternate circuit.

---

## 6. Input ownership and Fast copy contract

`EvaluateNew` must return an independently-owned ciphertext result.

Do not use full ciphertext `CopyNew()` on a Fast ciphertext with dormant high limbs.

The Mod1 evaluator may create a structural ciphertext at `LevelQ`, but initial cloning/copying must read only:

```text
q0/q1
metadata
```

Do not copy or synchronize q2...qL.

The input should remain usable after the call. Prefer not to mutate its q0/q1 coefficients or metadata.

Since the intended Fast C2S -> EvalMod handoff should already produce exactly the Mod1 starting Level, 007B may require:

```text
ct.Level() == Mod1Parameters.LevelQ
```

rather than implementing the Standard `ct.Level() > LevelQ` DropLevel convenience path.

If current C2S source evidence proves a broader Level range is required, support structural DropLevel only; never perform full-Q arithmetic to do so.

---

## 7. Reuse Fast polynomial evaluation directly

After normalization and cosine offset, evaluate:

```text
Mod1Parameters.Mod1Poly
```

through the accepted 007A Fast polynomial evaluator.

Do not create another PowerBasis or polynomial engine inside `mod1`.

Do not call Standard CKKS polynomial evaluation.

For the required `EvaluateNew` path, `scaling = 1`, so preserve Standard coefficient/scaling semantics exactly.

Avoid unnecessary per-call cloning of the Mod1 polynomial coefficients if source inspection confirms the Fast polynomial evaluator does not mutate the polynomial. If immutability is not sufficiently clear, a per-call clone is acceptable for correctness; report its allocation cost rather than introducing unsafe shared mutation.

---

## 8. Implement DoubleAngle with Fast primitives

After Fast polynomial evaluation, preserve the exact Standard loop:

```text
sqrt2pi = Sqrt2Pi

repeat DoubleAngle times:
    sqrt2pi = sqrt2pi * sqrt2pi
    res = MulRelin(res, res)
    res = Add(res, res)            // 2 * res^2
    res = Add(res, -sqrt2pi)       // 2 * res^2 - a^2
    res = Rescale(res)
```

Use only:

```text
Fast MulRelin
Fast Add
Fast Rescale
```

This is the homomorphic implementation of the amplitude-aware cosine double-angle recurrence.

Do not use Standard Relinearize, evaluation keys, GadgetProduct, QP, or full-Q arithmetic.

Only q0/q1 may influence the result.

---

## 9. Restore output Scale semantics exactly

The Standard Mod1 source ends with:

```go
res.Scale = ct.Scale
```

for the `EvaluateNew` path. Preserve the same metadata meaning in Fast Mod1.

Then add the Bootstrap-level Fast boundary:

```go
func (eval *FastEvaluator) EvalMod(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, error)
```

which delegates to Fast Mod1 and then, exactly like Standard Bootstrap:

```go
ctOut.Scale = eval.Parameters.BootstrappingParameters.DefaultScale()
```

Do not add an arithmetic multiply merely because comments describe this as “multiply back by q”; the current Standard implementation performs the final restoration through Scale metadata at this boundary.

---

## 10. Output contract

Fast Bootstrap `EvalMod` output must be suitable for the existing Fast SlotsToCoeffs path:

```text
Degree = 1
NTT = true
Montgomery = true
q0/q1 authoritative
q2...qL dormant
Scale = BootstrappingParameters.DefaultScale()
Level = the Standard-equivalent post-Mod1 level for the same parameters
LogDimensions / batching metadata preserved
```

The exact output Level must be derived from the current Mod1 polynomial depth + DoubleAngle rescale schedule, not hardcoded.

Do not require dormant q2...qL equality with Standard.

---

## 11. Correctness oracle strategy

Use a real numerical oracle for the complete Mod1 path.

### A. Standard Mod1 reference

Preferred primary reference:

1. Construct the same `mod1.Parameters` from the current literal.
2. Build a Standard CKKS evaluator + polynomial evaluator + Mod1 evaluator in test code only.
3. Create deterministic, valid full-Q zero-secret-style input:
   - encoded message in c0;
   - c1 = 0;
   - sufficient Level;
   - valid NTT/Montgomery representation.
4. Run Standard Mod1 on an authoritative full-Q copy.
5. Run Fast Mod1 on a copy whose q0/q1 are identical but q2...qL may be independently poisoned.
6. Decode both results through appropriate test/reference boundaries.
7. Compare decoded numerical values within a documented CKKS tolerance.

A Standard relinearization key may be generated in test code if required by the Standard evaluator. Do not import this requirement into production Fast code.

This comparison validates the complete current-source normalization + polynomial + DoubleAngle semantics without inventing a second Mod1 formula.

### B. Clear-function sanity oracle

Also include at least one small deterministic test checking the intended periodic behavior in the supported input domain, conceptually:

```text
x = integer + epsilon
Fast Mod1(x) ≈ epsilon
```

Account for the exact Standard normalization assumptions (`K`, QDiff/correction, Scale) rather than comparing against a naive `% 1` on an incorrectly-scaled ciphertext.

If a direct clear `%1` comparison would require reproducing hidden C2S correction factors, use a clear evaluation of the exact generated Mod1 polynomial + DoubleAngle circuit instead and document the boundary.

---

## 12. Poisoned dormant-residue test

Use two Fast inputs with identical:

```text
q0/q1
metadata
```

but different poison in:

```text
q2...qL
```

Run the complete Fast Mod1 / Bootstrap EvalMod path.

Verify:

```text
q0/q1 outputs identical
Level/Scale/representation identical
poison does not influence the answer
no deliberate q2...qL synchronization is required
```

This test should fail if any production path silently routes through Standard full-RNS arithmetic.

---

## 13. Required tests

Add at least:

### A. Constructor / parameters

Verify Fast Bootstrap constructs the complete Mod1 parameters and the Fast polynomial/Mod1 evaluators from the same literals used by Standard.

### B. Input validation

Reject:

```text
nil evaluator/input/metadata
wrong ring type
wrong degree
wrong Level for the bounded handoff
non-NTT
non-Montgomery
unsupported Mod1Type
non-nil Mod1InvPoly / Mod1InvDegree > 0
```

### C. normalization + offset

Use a small controlled configuration or test hook to ensure omission of the Scale normalization or cosine offset would make the test fail.

Do not merely check that the code runs.

### D. DoubleAngle

Use a configuration with nonzero DoubleAngle and verify the complete numerical result against Standard/reference.

The target/default-like case with `DoubleAngle = 3` must be exercised.

### E. full CosDiscrete Mod1 numerical oracle

Use a default-like:

```text
K = 16
Mod1Degree ≈ 30
DoubleAngle = 3
Mod1InvDegree = 0
```

and compare decoded Fast vs Standard/reference output.

### F. Bootstrap EvalMod wrapper

Verify final:

```text
Scale = BootstrappingParameters.DefaultScale()
Level matches Standard-equivalent schedule
NTT/Montgomery preserved
Degree = 1
```

### G. poisoned high limbs

As specified above.

### H. repeated evaluation

Reuse the same Fast Mod1/polynomial evaluator for multiple inputs and verify no workspace contamination of earlier returned results.

---

## 14. Performance benchmarks

Benchmark warm steady-state Fast EvalMod at:

```text
LogN = 13
LogN = 16
```

Use a default-like CosDiscrete configuration with:

```text
K = 16
Mod1Degree ~30
DoubleAngle = 3
Mod1InvDegree = 0
```

The benchmark must include:

```text
normalization/offset
Fast polynomial evaluation
all DoubleAngle rounds
final Scale restoration
```

Report:

```text
ns/op
B/op
allocs/op
```

Warm evaluator/polynomial workspace before timing.

A Standard full-Q Mod1 benchmark is encouraged if it can be set up fairly without distorting the task, but do not require a hard ratio for acceptance.

Interpret the result against the already-known 007A polynomial cost. If most EvalMod time is inherited from Fast polynomial evaluation, report that instead of inventing a premature micro-optimization task.

Do not open 007B-P solely because structural ciphertext allocation remains high unless profiling shows an accidental full-Q/dormant-limb regression.

---

## 15. Preserve existing tests

Run at least:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/polynomial
go test ./circuits/ckks/mod1
go test ./circuits/ckks/dft
go test ./circuits/ckks/bootstrapping
```

Do not regress:

```text
Fast polynomial evaluation
Fast Rescale
Fast DFT
Fast ScaleDown
Fast ModUp / Trace
Standard Mod1
Standard Bootstrap
```

---

## 16. Out of scope

Do not implement in 007B:

```text
full Fast Bootstrap orchestration
automatic CoeffsToSlots -> EvalMod -> SlotsToCoeffs wiring
SinContinuous
CosContinuous
ArcSine / Mod1InvPoly
EvalModAndScale for iterative bootstrapping
DenseToSparse / SparseToDense
QP arithmetic
GadgetProduct
Standard key switching
new ciphertext storage representation
noise-fidelity injection
```

Do not modify Standard Mod1 mathematics or Standard CKKS behavior.

---

## 17. Documentation

Update `docs/FAST_CKKS_SPEC.md` after implementation so status becomes conceptually:

```text
Fast polynomial execution: implemented
Fast EvalMod / Mod1: implemented for Stage-A CosDiscrete, no-inverse path
full Fast Bootstrap orchestration: not implemented
```

Clearly state the bounded unsupported cases rather than implying generic Mod1 support.

---

## Preferred files

Prefer a narrow change set around:

```text
circuits/ckks/mod1/fast.go
circuits/ckks/mod1/fast_test.go
circuits/ckks/mod1/fast_bench_test.go
circuits/ckks/bootstrapping/fast_scaledown.go   // constructor/fields only if needed
circuits/ckks/bootstrapping/fast_evalmod.go
circuits/ckks/bootstrapping/fast_evalmod_test.go
docs/FAST_CKKS_SPEC.md
```

Small test-only Standard-reference setup is allowed.

Avoid modifying:

```text
circuits/common/polynomial
circuits/ckks/polynomial/fast.go
schemes/schemes.go
core/rlwe
```

unless current source proves a small fix is genuinely required. If 007A must be changed, report the source-backed reason explicitly rather than silently broadening scope.

---

## Completion

Commit with:

```text
feat(ckks/fast): implement Fast EvalMod
```

Push to:

```text
origin/fast-ckks
```

Final report must include:

1. Fast Mod1 API/location;
2. supported Mod1Type/inverse/scaling scope;
3. complete Mod1 parameter construction;
4. input Level/domain contract;
5. q0/q1-only input copy strategy;
6. exact Scale normalization behavior;
7. cosine offset/change-of-variable behavior;
8. polynomial evaluator reuse;
9. DoubleAngle implementation and count;
10. output Scale restoration;
11. Bootstrap `FastEvaluator.EvalMod` wrapper;
12. output Level/Scale/domain;
13. Standard/reference numerical comparison and tolerance;
14. periodic/clear sanity result;
15. poisoned q2...qL result;
16. repeated-evaluation result;
17. LogN13 benchmark;
18. LogN16 benchmark;
19. B/op and allocs/op;
20. tests run;
21. files changed;
22. commit hash;
23. push result.
