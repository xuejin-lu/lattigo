# Task 009 — End-to-End Fast Bootstrap

## Goal

Wire the already-implemented Fast stages into the first complete Stage-A Fast CKKS Bootstrap path.

Base reviewed commit:

```text
4fdb44a23a185a181d132256e419d54a28ec2c6e
feat(ckks/fast): add Fast bootstrap packing boundary
```

The target orchestration is:

```text
Fast Bootstrap / BootstrapMany
    ↓
PackAndSwitchN1ToN2
    ↓
ScaleDown
    ↓
ModUp + Trace
    ↓
Fast CoeffsToSlots
    ↓
Fast EvalMod
    ↓
Fast SlotsToCoeffs
    ↓
UnpackAndSwitchN2ToN1
    ↓
restore public residual Scale/domain
```

Task 009 is primarily an integration task. Do not redesign already-accepted Fast arithmetic, DFT, polynomial, EvalMod, packing, or ring-degree algorithms unless integration exposes a concrete source-backed defect.

Task 010 remains the dedicated final profiling/cleanup stage.

---

## 1. Read current source first

Before editing, read at least:

```text
circuits/ckks/bootstrapping/evaluator.go
circuits/ckks/bootstrapping/parameters.go
circuits/ckks/bootstrapping/bootstrapper.go
circuits/ckks/bootstrapping/fast_scaledown.go
circuits/ckks/bootstrapping/fast_modup.go
circuits/ckks/bootstrapping/fast_evalmod.go
circuits/ckks/bootstrapping/fast_packing.go
circuits/ckks/dft/fast.go
circuits/ckks/mod1/fast.go
schemes/ckks/fast/ring_degree.go
schemes/ckks/fast/partial_ntt.go
docs/FAST_CKKS_SPEC.md
```

Preserve the current Standard Bootstrap order and parameter semantics.

---

## 2. Important integration facts discovered before 009

### A. Fast DFT execution exists but Bootstrap does not yet own DFT matrices

`circuits/ckks/dft/fast.go` already implements Fast:

```text
CoeffsToSlots
SlotsToCoeffs
```

but `bootstrapping.FastEvaluator` currently does not own:

```text
*dft.FastEvaluator
C2SDFTMatrix
S2CDFTMatrix
```

and its constructor does not perform the Standard Bootstrap DFT scaling preparation.

### B. Standard Bootstrap initialization modifies DFT scaling before matrix generation

The current Standard `Evaluator.initialize` computes the Mod1-dependent factors conceptually as:

```text
qDiv = ScalingFactor / rounded(Q0)
qDiv = min(qDiv, 1)

C2S scaling *= qDiv / (K * QDiff)
S2C scaling *= DefaultScale / (ScalingFactor / MessageRatio)
```

and only then calls `dft.NewMatrixFromLiteral`.

These factors are part of the Bootstrap mathematics. Fast must not generate raw DFT matrices from the unadjusted literals.

### C. Standard Bootstrap accepts Level-0 input

Current Standard Bootstrap tests explicitly encrypt a Level-0 plaintext/ciphertext and call `Bootstrap`.

However the current Task-008 Fast packing boundary and current Fast ring-degree primitive require q0/q1 (`Level >= 1`).

This is a real integration gap. 009 must support the current Stage-A Level-0 Bootstrap case rather than silently requiring Level >= 1.

### D. Fast core becomes Montgomery after ModUp

Current Fast ModUp converts maintained residues to Montgomery and all Fast DFT/EvalMod execution is intentionally NTT + Montgomery.

A public Bootstrap result must not accidentally remain in an internal representation that ordinary decrypt/decode code interprets incorrectly.

### E. Dormant high limbs are not decryptable through the ordinary full-RNS decoder

A public result at Level > 1 whose q2...qL values are dormant is not a standard fully-authoritative ciphertext.

For Stage A, do not pretend otherwise.

---

## 3. Bounded Stage-A public Bootstrap scope

Required support:

```text
ring.Standard
CircuitOrder = ModUpThenEncode
IterationsParameters == nil
EphemeralSecretWeight == 0
Mod1Type = CosDiscrete
Mod1InvPoly == nil / Mod1InvDegree = 0
single non-iterative EvalMod path
Degree = 1
NTT input
ordinary/non-Montgomery public input
N2 == N1 or N2 == 2*N1
input Level >= 0
```

Required current public-output bound:

```text
ResidualParameters.MaxLevel() <= 1
```

Reason: current Fast public result maintains q0 at Level 0 and q0/q1 at Level 1. Returning a logical Level > 1 result with dormant q2+ as if it were an ordinary fully-authoritative ciphertext would break ordinary decrypt/decode semantics.

If the actual current engineering target has `ResidualParameters.MaxLevel() <= 1`, use it directly.

If a test/configuration requires residual MaxLevel > 1, do **not** synchronize dormant high limbs merely to make the test pass. Return an explicit unsupported-Stage-A error unless a deliberate final materialization boundary is proven necessary and is separately justified. Do not broaden 009 into full-RNS restoration.

Out of scope:

```text
ConjugateInvariant
SinContinuous
CosContinuous
ArcSine / Mod1InvPoly
iterative/META bootstrap
PREC128 iterative correction path
DenseToSparse / SparseToDense
arbitrary N2/N1 ratio
full dormant-limb materialization
noise-fidelity work
```

---

## 4. Preserve minimal FastEvaluator construction

Existing tests instantiate `NewFastEvaluator` with deliberately partial parameter objects for ScaleDown, ModUp, EvalMod, and packing unit tests.

Do not break those narrow constructions merely because end-to-end Bootstrap now needs DFT matrices.

Preferred design:

```text
NewFastEvaluator
    -> constructs existing cheap Fast primitives
    -> full Bootstrap DFT/matrix initialization remains lazy

Bootstrap / BootstrapMany / bootstrapCore
    -> ensureFastBootstrapCircuit()
```

Equivalent eager construction is acceptable only if all existing minimal constructors continue to work without inventing fake DFT/Mod1 parameters.

---

## 5. Add complete Fast Bootstrap circuit data

Extend `FastEvaluator` with the bounded full-circuit state, conceptually:

```go
DFTEvaluator  *dft.FastEvaluator
C2SDFTMatrix  dft.Matrix
S2CDFTMatrix  dft.Matrix

fastBootstrapInitialized bool
fastBootstrapErr error
```

Naming may differ.

### Required initialization parity

Fast must use exactly the current Standard Bootstrap parameter preparation for:

```text
Mod1Parameters
C2SScaling
S2CScaling
C2S matrix generation
S2C matrix generation
```

Preferred implementation:

Extract the existing pure parameter/matrix-building portion of Standard `initialize` into a package-local helper used by both Standard and Fast evaluators.

For example conceptually:

```go
func buildBootstrapCircuitData(params Parameters) (
    adjusted Parameters,
    mod1Params mod1.Parameters,
    c2s dft.Matrix,
    s2c dft.Matrix,
    err error,
)
```

The exact API is flexible.

This helper must perform planning/encoding only. It must not construct or invoke a Standard CKKS evaluator, key switching, QP, or GadgetProduct.

If a shared helper would require a risky refactor, a Fast-local source-faithful mirror is acceptable, but then add an explicit parity test against current Standard initialization for the scaling values and matrix schedule.

Do not silently change Standard numerical behavior.

---

## 6. Level-0 maintained-limb support at the packing boundary

009 must close the integration gap discovered after Task 008.

General Fast maintained-limb rule:

```text
Level 0  -> q0 authoritative
Level >=1 -> q0/q1 authoritative
```

### A. Fast packing/unpacking

Generalize Task-008 helpers from fixed `limb < 2` to the maintained limb count:

```go
maintained := min(Level+1, 2)
```

At Level 0:

```text
copy q0 only
pack monomial multiply/add on q0 only
unpack monomial multiply on q0 only
```

At Level >= 1, current q0/q1 behavior remains unchanged.

Monomial tables for a ring with MaxLevel 0 must be allowed to contain q0 only.

Existing Level>=1 tests must continue to pass unchanged.

### B. Fast ring-degree conversion

Generalize `FastN1ToN2` / `FastN2ToN1` to support Level 0.

For NTT Level-0 conversion, do not route through the current two-limb-only `FastPartialNTT/FastPartialINTT` validator. Either:

```text
- add a maintained-limb transform helper that handles 1 or 2 limbs; or
- use direct SubRing[0] INTT/NTT locally for the Level-0 path.
```

Do not weaken the existing q0/q1 partial-transform contract if a local ring-degree solution is cleaner.

Add Level-0 coefficient-domain and NTT-domain ring-degree tests.

No evaluation key may be introduced.

---

## 7. Public input contract

`Bootstrap` / `BootstrapMany` must reject early with clear errors when the supported Stage-A contract is violated.

At minimum validate:

```text
non-nil evaluator/input/metadata
non-empty BootstrapMany slice
Standard ring
Degree = 1
IsNTT = true
IsMontgomery = false
input N matches ResidualParameters.N()
input Level within ResidualParameters
same Level across BootstrapMany inputs
same LogSlots / batching representation where packing requires it
CircuitOrder = ModUpThenEncode
IterationsParameters == nil
EphemeralSecretWeight == 0
supported Mod1 configuration
ResidualParameters.MaxLevel() <= 1
N2 == N1 or N2 == 2*N1
```

Why reject Montgomery input at the public boundary:

```text
ScaleDown can preserve it,
but current Fast ModUp basis explicitly requires ordinary/non-Montgomery input.
```

Fail early instead of failing halfway through Bootstrap.

---

## 8. Implement the Fast Bootstrap core

Prefer a private method conceptually:

```go
func (eval *FastEvaluator) bootstrapCore(
    ctIn *rlwe.Ciphertext,
) (ctOut *rlwe.Ciphertext, errScale *rlwe.Scale, err error)
```

Mirror current Standard `bootstrap` stage order exactly:

```text
1. ScaleDown
2. ModUp
3. CoeffsToSlots
4. EvalMod(real)
5. EvalMod(imag), if present
6. SlotsToCoeffs
```

Pseudo-structure:

```go
ct, errScale, err := eval.ScaleDown(ctIn)
ct, err = eval.ModUp(ct)
ctReal, ctImag, err := eval.DFTEvaluator.CoeffsToSlotsNew(ct, eval.C2SDFTMatrix)
ctReal, err = eval.EvalMod(ctReal)
if ctImag != nil {
    ctImag, err = eval.EvalMod(ctImag)
}
ctOut, err = eval.DFTEvaluator.SlotsToCoeffsNew(ctReal, ctImag, eval.S2CDFTMatrix)
```

Do not add a Standard evaluator fallback.

Do not use evaluation keys in the Fast core.

The `errScale` value may be returned for source parity but the non-iterative Stage-A public path does not need the iterative correction branch.

---

## 9. Stage transition contracts

Add assertions/tests for the real integration boundaries rather than only testing final output.

For a valid parameter set, verify conceptually:

### After ScaleDown

```text
Level = 0
NTT = true
Montgomery = false
Scale approximately q0 / MessageRatio
```

### After ModUp

```text
Level = BootstrappingParameters.MaxLevel()
NTT = true
Montgomery = true
q0/q1 authoritative
q2+ dormant
```

### After CoeffsToSlots

```text
Level = Mod1Parameters.LevelQ
NTT/Montgomery preserved
ctImag presence follows the DFT format / LogSlots case
```

### After EvalMod

```text
Level = SlotsToCoeffsParameters.LevelQ
Scale = BootstrappingParameters.DefaultScale()
NTT/Montgomery preserved
```

Use the exact source-derived level schedule; do not hardcode these equalities if current parameter methods expose a more precise expression.

### After SlotsToCoeffs

```text
Level = ResidualParameters.MaxLevel()
NTT = true
Montgomery = true before the public finalization boundary
```

If any stage's actual source contract differs, preserve the source and document the difference rather than forcing these comments literally.

---

## 10. Public output finalization

After the Fast core and after unpack/ring-switch restoration, return an ordinary public ciphertext compatible with the current Stage-A residual output.

Required:

```text
N = ResidualParameters.N()
Degree = 1
Level = ResidualParameters.MaxLevel()   // currently bounded to 0 or 1
Scale = ResidualParameters.DefaultScale()
IsNTT = true
IsMontgomery = false
original LogSlots restored
batch/bit-reverse metadata preserved consistently
```

### Montgomery exit

Before public return, convert every active maintained limb:

```text
Level 0  -> q0 IMForm
Level 1  -> q0/q1 IMForm
```

for c0 and c1, then set:

```go
IsMontgomery = false
```

Do not IMForm dormant q2+.

This is a deliberate public representation boundary, not a general full-Q synchronization step.

---

## 11. Implement Bootstrap and BootstrapMany

Add:

```go
func (eval *FastEvaluator) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error)

func (eval *FastEvaluator) BootstrapMany(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error)
```

`BootstrapMany` must mirror Standard outer orchestration:

```text
PackAndSwitchN1ToN2
→ bootstrapCore on each packed ciphertext
→ UnpackAndSwitchN2ToN1
→ residual DefaultScale
→ public Montgomery exit
```

`Bootstrap` may wrap `BootstrapMany` for one ciphertext.

Do not duplicate the core circuit between the two APIs.

### Interface methods

Implement the remaining `Bootstrapper` surface if absent:

```go
Depth() int
MinimumInputLevel() int
OutputLevel() int
```

The Fast evaluator should satisfy:

```go
var _ Bootstrapper = (*FastEvaluator)(nil)
```

Use values that describe the actual supported Fast path. In particular, Level-0 input is supported in Stage A, so `MinimumInputLevel()` should not falsely reject the supported public case.

---

## 12. Input ownership / mutation

Do not require callers to sacrifice their source ciphertext merely because the Fast stages mutate internal buffers in place.

At the public `Bootstrap` / `BootstrapMany` boundary, ensure the working ciphertext storage is independent from caller-owned q0/q1 data.

Do not solve this with full `CopyNew()` on dormant high limbs.

Reuse the maintained-copy pattern:

```text
allocate required structural ciphertext
copy metadata
copy q0 at Level0
copy q0/q1 at Level>=1
```

Packing already owns/copies outputs; when a no-op packing case would otherwise alias the caller, preserve independent working ownership explicitly.

Add a test that the source ciphertext remains unchanged after Bootstrap.

---

## 13. End-to-end numerical correctness

This is the most important acceptance criterion.

### A. Clear message preservation oracle

Use deterministic nonzero CKKS slot values.

Construct a valid current Stage-A zero-secret-style input:

```text
c0 = encoded plaintext
c1 = 0
```

with the public input metadata/domain expected by Fast Bootstrap.

Run:

```text
Fast Bootstrap
→ ordinary public output
→ ordinary decrypt/decode or direct plaintext decode under zero-secret semantics
```

Compare decoded output against the original clear values within a documented CKKS tolerance.

Do not use all-zero input as the main oracle.

### B. Standard Bootstrap numerical reference

For a small/insecure test parameter set, also run current Standard Bootstrap on an equivalent valid zero-secret-style/full-Q input and compare decoded numerical values:

```text
Fast decoded output
vs
Standard decoded output
```

within a documented tolerance.

Standard test-only evaluation keys are allowed.

Production Fast code must not depend on them.

This reference should verify the complete source-level combination of:

```text
C2S scaling
EvalMod normalization
S2C scaling
final output scale
```

### C. Do not use raw ciphertext equality as the end-to-end oracle

Fast intentionally does not reproduce Standard key-switch/noise representation.

Use decoded semantics and maintained-residue boundary tests where appropriate.

---

## 14. Required end-to-end cases

### Case 1 — Same ring, single ciphertext, Level 0

This is mandatory because current Standard Bootstrap explicitly supports Level-0 input.

Verify:

```text
Bootstrap succeeds
message preserved numerically
output Level = residual MaxLevel
output Scale = residual DefaultScale
output IsMontgomery = false
ordinary decrypt/decode succeeds
```

### Case 2 — Same ring, Level 1 input

Exercise the q0/q1 input path and ScaleDown transition.

### Case 3 — Full-slot C2S branch

Use `LogSlots == LogMaxSlots` so `CoeffsToSlots` produces both real and imaginary ciphertexts.

Verify both receive EvalMod and recombine correctly in SlotsToCoeffs.

### Case 4 — Sparse/repacked C2S branch

Use `LogSlots < LogMaxSlots` / `RepackImagAsReal` so `ctImag == nil` after C2S.

Verify the one-ciphertext EvalMod path.

### Case 5 — BootstrapMany odd count

Use at least 3 inputs and verify order, output count, restored LogSlots, and independent ownership.

### Case 6 — Factor-two N1→N2→N1

Use `N2 = 2*N1` and include Level-0 coverage after the ring-degree generalization.

Verify decoded semantics after full:

```text
pack N1
→ ring switch
→ Fast bootstrap core
→ ring switch back
→ unpack N1
```

If the small factor-two numerical case is excessively expensive, a stage-correctness integration test plus a same-ring numerical full Bootstrap is acceptable only if the factor-two outer boundary remains covered by exact q0/q1 reference tests from Task 008 plus new Level-0 ring-degree tests. Report the limitation explicitly.

---

## 15. Dormant-residue independence across the full core

Where a high-level Bootstrap public test cannot contain q2+ because `ResidualParameters.MaxLevel() <= 1`, add an internal core test at a higher Bootstrap input level:

```text
same q0/q1
same metadata
different q2...qL poison
```

Run through as much of:

```text
ScaleDown
→ ModUp
→ C2S
→ EvalMod
→ S2C
```

as the source-valid level schedule permits.

Verify maintained q0/q1 output is identical.

Do not manufacture an invalid public residual configuration merely to create poison rows.

Existing per-stage poison tests remain necessary but are not a substitute for at least one integrated poison test if practical.

---

## 16. DFT initialization parity tests

If Standard and Fast share the same extracted pure initialization helper, test that both constructors receive the same:

```text
adjusted C2S scaling
adjusted S2C scaling
Mod1 parameters
C2S LevelQ / level schedule
S2C LevelQ / level schedule
```

If the logic is mirrored instead of shared, this parity test is mandatory.

At least one numerical end-to-end test must fail if either C2S or S2C scaling is intentionally omitted.

---

## 17. Validation tests

Explicitly reject:

```text
nil ciphertext
empty BootstrapMany
wrong ring type
Degree != 1
non-NTT public input
Montgomery public input
unsupported CircuitOrder
IterationsParameters != nil
EphemeralSecretWeight != 0
unsupported Mod1 type/inverse
Residual MaxLevel > 1 without materialization support
N2 < N1
N2 > 2*N1
mismatched input levels in BootstrapMany
mismatched LogSlots where packing requires consistency
```

Do not silently route any unsupported case to Standard.

---

## 18. Performance benchmark

Add warm end-to-end Fast Bootstrap benchmarks.

Required if practical:

```text
LogN = 13
LogN = 16
```

Use the actual current engineering parameter profile if available in the repository/test helpers. Otherwise use a source-consistent default-like parameter set with:

```text
CosDiscrete
Mod1 degree about 30
DoubleAngle = 3
no inverse
same-ring single ciphertext
```

Benchmark must include the public path:

```text
packing boundary
ScaleDown
ModUp
C2S
EvalMod
S2C
unpacking boundary
public Montgomery exit
```

Warm:

```text
DFT matrices
packing monomial tables
polynomial workspace
other evaluator-owned reusable state
```

before `ResetTimer()`.

Report:

```text
ns/op
B/op
allocs/op
input Level
output Level
LogN
LogSlots
```

A Standard end-to-end benchmark is encouraged if setup is fair and runtime is practical, but no hard speed ratio is required for 009.

Do not create 009-P solely because allocation is still structurally high. Task 010 is the profiling/cleanup stage unless 009 reveals a correctness bug or accidental Standard/full-Q fallback.

---

## 19. No hidden Standard execution

Production Fast Bootstrap must not call:

```text
ckks.NewEvaluator for execution
Standard dft.Evaluator
Standard mod1.Evaluator
ApplyEvaluationKey
GadgetProduct
QP arithmetic
Standard key switching
Standard relinearization
full-Q arithmetic on dormant rows
```

Using shared **pure planning/matrix-generation** code is allowed and preferred.

Test-only Standard reference execution is allowed.

---

## 20. Regression commands

Run at least:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/polynomial
go test ./circuits/ckks/mod1
go test ./circuits/ckks/dft
go test ./circuits/ckks/bootstrapping
```

Also run targeted end-to-end Fast Bootstrap tests and the new benchmarks.

No regressions to Tasks 003–008.

---

## 21. Documentation

Update `docs/FAST_CKKS_SPEC.md` after implementation.

The status should become conceptually:

```text
Fast ScaleDown              implemented
Fast ModUp + Trace          implemented
Fast C2S                    implemented
Fast EvalMod                implemented
Fast S2C                    implemented
Fast Pack/Unpack/ring       implemented
Fast Bootstrap orchestration implemented for bounded Stage-A Standard-ring path
```

Document the exact current public restrictions:

```text
CosDiscrete/no inverse
non-iterative
no ephemeral switching
N2 == N1 or 2*N1
Residual MaxLevel <= 1 unless a deliberate materialization boundary is later added
```

Do not imply generic Standard-Bootstrap feature parity.

---

## Preferred files

Prefer changes around:

```text
circuits/ckks/bootstrapping/fast_bootstrap.go
circuits/ckks/bootstrapping/fast_bootstrap_test.go
circuits/ckks/bootstrapping/fast_bootstrap_bench_test.go
circuits/ckks/bootstrapping/fast_scaledown.go
circuits/ckks/bootstrapping/fast_packing.go
schemes/ckks/fast/ring_degree.go
schemes/ckks/fast/ring_degree_test.go
docs/FAST_CKKS_SPEC.md
```

A small pure initialization helper extracted from Standard `evaluator.go` is allowed if it is behavior-preserving and shared by Standard/Fast.

Avoid changing already-accepted arithmetic/DFT/Mod1 code unless a failing integration oracle proves it necessary.

---

## Completion

Commit with:

```text
feat(ckks/fast): wire end-to-end Fast Bootstrap
```

Push to:

```text
origin/fast-ckks
```

Final report must include:

1. `Bootstrap` / `BootstrapMany` APIs and interface status;
2. supported Stage-A configuration;
3. DFT matrix/scaling initialization strategy;
4. confirmation Standard initialization parity;
5. Level-0 packing support change;
6. Level-0 ring-degree support change;
7. complete stage order;
8. stage Level/Scale/domain transitions;
9. full-slot real+imag path result;
10. sparse/repacked path result;
11. same-ring Level-0 numerical result and tolerance;
12. Level-1 numerical result;
13. Standard-vs-Fast decoded reference result;
14. BootstrapMany odd-count result;
15. factor-two boundary result;
16. input ownership result;
17. dormant-poison integrated result;
18. final public output Level/Scale/NTT/Montgomery state;
19. proof ordinary decrypt/decode works on public output;
20. confirmation no Standard/QP/Gadget/key-switch production fallback;
21. LogN13 benchmark;
22. LogN16 benchmark if practical;
23. B/op and allocs/op;
24. remaining limitations;
25. regression tests;
26. files changed;
27. commit hash;
28. push result.
