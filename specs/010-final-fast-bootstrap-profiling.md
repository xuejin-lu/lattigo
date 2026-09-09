# Task 010 — Final Fast Bootstrap Profiling and Targeted Optimization

## Goal

Profile the completed Stage-A Fast Bootstrap end to end, identify the real CPU/allocation bottlenecks, and apply only high-value optimizations that preserve the accepted Fast architecture and numerical behavior.

Reviewed base:

```text
29a48a54ddd89c6ec3a2ea6cabe0cfa09df06b18
test(ckks/fast): close end-to-end Bootstrap review gaps
```

Task 009/009-P is accepted. Do not reopen Bootstrap functionality unless profiling exposes a real correctness/performance-path defect.

This is the final profiling task, not another feature expansion.

---

## 1. Preserve the accepted functional contract

The following path is already accepted:

```text
BootstrapMany
→ Fast pack / optional N1→N2
→ ScaleDown
→ ModUp + Trace
→ Fast CoeffsToSlots
→ Fast EvalMod
→ Fast SlotsToCoeffs
→ optional N2→N1 / unpack
→ public IMForm finalization
```

Preserve:

```text
ring.Standard
CosDiscrete
no Mod1 inverse
non-iterative Stage-A path
N2=N1 or N2=2*N1
Level0/Level1 public input
Residual MaxLevel <= 1 public output
q0 at Level0; q0/q1 at Level>=1
no Standard execution fallback
no QP/Gadget/key switching
```

All existing numerical oracles from 003–009-P must continue to pass.

---

## 2. Current measured baseline

009-P reported on Apple M4 with:

```text
LogSlots = 4
input Level = 0
output Level = 1
K = 16
Mod1 degree = 30
DoubleAngle = 3
Mod1 inverse = 0
-benchtime=1x
```

Results:

```text
LogN13
Fast      ~11.38 ms,  33.6 MB/op,  3947 allocs/op
Standard ~187.24 ms,  64.5 MB/op, 21988 allocs/op

LogN16
Fast     ~123.96 ms, 279.7 MB/op,  4301 allocs/op
Standard   ~2.200 s, 510.5 MB/op, 24436 allocs/op
```

These establish a strong performance trend but are not stable statistics because `-benchtime=1x` was used.

Task 010 must produce a more stable final baseline and post-optimization measurement.

---

## 3. Phase A — stable end-to-end baseline

Use the existing fair Fast/Standard benchmark and keep identical:

```text
parameters
message
input Level
LogSlots
K
Mod1 degree
DoubleAngle
public API path
CopyNew treatment inside timed loops
```

Run at least:

```bash
go test ./circuits/ckks/bootstrapping \
  -run '^$' \
  -bench 'BenchmarkFastBootstrapEndToEnd' \
  -benchmem \
  -benchtime=3x \
  -count=3
```

If LogN16 Standard makes this disproportionately slow, `-benchtime=2x -count=3` is acceptable, but use the same setting for Fast and Standard within that LogN.

Record per run:

```text
ns/op
B/op
allocs/op
```

Report median or clearly labeled representative values; do not present a single noisy run as exact performance.

---

## 4. Phase B — collect CPU and allocation profiles before editing

Profile **Fast only**, preferably LogN13 first because it is cheap enough for repeated samples.

Use Go benchmark profiling, for example:

```bash
go test ./circuits/ckks/bootstrapping \
  -run '^$' \
  -bench 'BenchmarkFastBootstrapEndToEnd/LogN=13/Fast$' \
  -benchmem \
  -benchtime=20x \
  -cpuprofile /tmp/fast-bootstrap-cpu.out \
  -memprofile /tmp/fast-bootstrap-mem.out

go tool pprof -top /tmp/fast-bootstrap-cpu.out
go tool pprof -top -alloc_space /tmp/fast-bootstrap-mem.out
go tool pprof -top -alloc_objects /tmp/fast-bootstrap-mem.out
```

Equivalent commands are acceptable.

For LogN16, use a smaller iteration count if needed.

Do not commit binary pprof files.

Create/update a text profiling record, preferably:

```text
docs/FAST_CKKS_PROFILE.md
```

Record the top CPU and allocation entries with enough context to identify the source path.

---

## 5. High-priority hypothesis to verify — dormant backing-array allocation

Fast maintains only:

```text
Level0   -> q0
Level>=1 -> q0/q1
```

but several accepted helpers still allocate ciphertexts with:

```go
ckks.NewCiphertext(params, degree, logicalLevel)
```

Standard `ring.NewPoly`, `ring.Poly.Resize`, `rlwe.NewElement`, and `rlwe.Element.Resize` allocate an `N`-length backing slice for every logical residue row.

Therefore a logical Level such as 10 can allocate q0...q10 even though Fast reads/writes only q0/q1.

Known candidate call sites include, but are not limited to:

```text
circuits/ckks/dft/fast.go
  CoeffsToSlotsNew
  SlotsToCoeffsNew

circuits/ckks/mod1/fast.go
  cloneFastCiphertext

circuits/ckks/polynomial/fast.go
  workspace initial allocation
  cloneMaintainedResult

circuits/ckks/bootstrapping/fast_packing.go
  copyFastPackingCiphertext
  ring-switch outputs

circuits/ckks/bootstrapping/fast_modup.go
  restoreFastModUpLevel

schemes/ckks/fast/evaluator_ntt.go
  *New methods / Resize paths
```

This is a hypothesis, not permission to rewrite all these files.

Only optimize it if `alloc_space` proves it is material.

---

## 6. Preferred optimization if dormant backing arrays are confirmed

If profiling confirms full logical-level ciphertext allocation is a major source of B/op, introduce a **localized Fast compact structural allocation strategy**.

Conceptual representation:

```text
logical Level = L
Coeffs slice length = L+1        // Level() remains L
q0 backing = length N
q1 backing = length N if L>=1
q2...qL backing = nil / absent data storage but structurally represented
```

Important:

```text
len(Coeffs) still encodes logical Level
only maintained rows own N-sized backing arrays
```

Do not change global Standard constructors such as:

```text
ring.NewPoly
rlwe.NewElement
ckks.NewCiphertext
```

Their behavior belongs to Standard Lattigo.

Prefer Fast-local helpers such as conceptually:

```go
newCompactFastCiphertext(...)
ensureFastMaintainedStorage(...)
resizeFastLogicalLevel(...)
```

Naming and location are flexible.

### Critical safety rule

Generic `ring.Poly.Resize` / `rlwe.Element.Resize` will allocate full backing rows when growing a level.

A compact Fast object must not accidentally pass through a grow operation that materializes q2...qL.

Audit every optimized path that changes Level.

Shrinking logical level by slicing headers is fine if ownership/later reuse remains correct.

Growing Level at deliberate ModUp must use the Fast compact logical-level helper rather than Standard grow behavior.

---

## 7. Compact-storage correctness contract

If compact internal storage is introduced, add explicit tests proving:

```text
Level() remains source-correct
N() remains source-correct
Degree() remains source-correct
q0/q1 arithmetic is unchanged
q2...qL are never read
nil dormant rows do not panic Fast execution
public Level0/1 output is fully materialized for its active maintained rows
ordinary zero-secret Decryptor + Decode still works
```

At a high internal level, assert structurally:

```text
len(poly.Coeffs) == Level+1
len(poly.Coeffs[0]) == N
len(poly.Coeffs[1]) == N   // when Level>=1
len(poly.Coeffs[i]) == 0   // i>=2, if compact representation is used
```

Do not serialize or pass compact internal ciphertexts to Standard/full-RNS APIs.

---

## 8. Other optimizations allowed only when profile-backed

After the first high-value allocation fix, reprofile before making another optimization.

Allowed candidates if they appear prominently in CPU/alloc profiles:

```text
- evaluator-owned reusable output ciphertexts/scratch
- removing avoidable maintained-result clones at private internal boundaries
- packing/unpacking temporary slice reuse
- ring-degree q0/q1 scratch reuse
- avoiding repeated big.Float/big.Int setup in hot loops
- avoiding repeated validation/planning work that is invariant per evaluator
```

Do not optimize low-ranked functions merely because they look inelegant.

Do not remove public ownership guarantees to save one allocation.

Do not introduce unsafe shared buffers across returned public ciphertexts.

---

## 9. CPU optimization discipline

For CPU time, identify top **self/cumulative** hot functions from pprof.

Classify them as:

```text
mathematical work that must remain
Fast-specific avoidable overhead
allocation/GC overhead
benchmark/setup overhead
```

Examples of required mathematical work that should not be removed solely for speed:

```text
q0/q1 NTT/INTT at real representation boundaries
Fast DFT linear transforms
polynomial multiplications/rescales
EvalMod DoubleAngle operations
public IMForm exit
```

Do not alter mathematics or reduce polynomial degree/K/DoubleAngle to improve benchmark numbers.

---

## 10. Stage attribution

Produce a high-level attribution of end-to-end Fast time/allocation across:

```text
packing / ring switch
ScaleDown
ModUp + Trace
CoeffsToSlots
EvalMod
SlotsToCoeffs
unpack / public finalization
```

Preferred evidence order:

1. pprof call stacks/top entries;
2. existing stage microbenchmarks under matching or nearby profiles;
3. new targeted stage benchmark only if needed to disambiguate.

Do not build a large instrumentation framework just for this task.

A percentage estimate may be approximate, but label estimates as estimates.

---

## 11. Reprofile after every material optimization

For each accepted optimization, capture:

```text
before ns/op, B/op, allocs/op
after  ns/op, B/op, allocs/op
percentage change
```

Use the same benchmark command/settings.

If an optimization reduces B/op strongly but makes CPU materially worse, do not automatically keep it. Judge against project priority:

```text
1. execution throughput
2. functional correctness
3. experiment flexibility
4. memory/allocation reduction
```

Memory work is valuable because it can improve GC and repeated 10k-run throughput, but speed remains primary.

---

## 12. Throughput-oriented final benchmark

The project goal is thousands of repeated inference/Bootstrap operations, not one cold call.

Add a Fast-only repeated benchmark or reuse the existing benchmark with enough iterations to measure warm steady state.

At minimum report stable Fast measurements for:

```text
LogN13
LogN16
```

Prefer:

```text
-benchtime=10x or higher for Fast
-count=3
```

when practical.

Do not include evaluator construction, DFT matrix generation, monomial-table initialization, or key generation inside timed steady-state loops.

---

## 13. Standard comparison

Keep Standard as the comparison baseline, but do not optimize Standard in this task.

Final report should include comparable Fast/Standard end-to-end numbers for LogN13 and LogN16.

Use wording such as:

```text
"On this Apple M4 benchmark/profile..."
```

Do not claim a universal speedup.

Do not compare Fast stage timing against a Standard benchmark that performs different mathematical work.

---

## 14. Functional regressions

After optimization, run at least:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/polynomial
go test ./circuits/ckks/mod1
go test ./circuits/ckks/dft
go test ./circuits/ckks/bootstrapping
```

The following 009-P oracles are mandatory survivors:

```text
full-slot real+imag numerical Bootstrap
Level-1 numerical Bootstrap
BootstrapMany distinct-message order
zero-secret ordinary Decryptor + Decode
Standard/Fast DFT planning parity
Standard-vs-Fast decoded reference
integrated dormant-poison independence
```

No tolerance widening solely to pass an optimization.

---

## 15. Architecture regression scan

Before completion, search changed production paths and verify no new production use of:

```text
Standard ckks.Evaluator execution
Standard dft.Evaluator execution
Standard mod1.Evaluator execution
ApplyEvaluationKey
GadgetProduct
QP arithmetic
full-Q arithmetic over dormant rows
per-coefficient arbitrary big.Int CRT in the hot path
```

Pure planning/matrix generation remains allowed.

---

## 16. Profiling document

Create/update:

```text
docs/FAST_CKKS_PROFILE.md
```

Keep it concise and evidence-based. Include:

```text
hardware / Go version if available
benchmark parameter profile
baseline table
CPU top hotspots
alloc_space top hotspots
selected optimization(s) and rationale
before/after table
remaining bottlenecks
known Stage-A restrictions
```

Do not paste full pprof dumps.

Update `docs/FAST_CKKS_SPEC.md` only if the implementation architecture actually changes, e.g. compact dormant-row storage becomes part of current Fast behavior.

---

## 17. Stop conditions

Do **not** keep optimizing indefinitely.

Stop after:

```text
- top allocation/CPU bottleneck(s) have been addressed or proven mathematically necessary;
- a second profile shows no obvious low-risk, high-value Fast-specific waste;
- functional tests pass;
- final stable benchmarks are recorded.
```

Do not create speculative micro-optimizations for tiny percentages.

If the top remaining cost is the actual q0/q1 DFT/EvalMod mathematics, report that clearly and stop.

---

## 18. Out of scope

Do not add:

```text
new Bootstrap features
new Mod1 modes
iterative Bootstrap
higher residual public levels through full-RNS materialization
security restoration
noise-fidelity calibration
application/CNN changes
parallel Bootstrap execution
GPU/FPGA acceleration
unsafe concurrent reuse of evaluator scratch
changes to Standard global storage semantics
```

Parallelism and hardware acceleration are separate future work, not Task 010.

---

## Preferred files

Profiling evidence:

```text
docs/FAST_CKKS_PROFILE.md
```

Potential production files only if profile-backed:

```text
schemes/ckks/fast/...
circuits/ckks/dft/fast.go
circuits/ckks/polynomial/fast.go
circuits/ckks/mod1/fast.go
circuits/ckks/bootstrapping/fast_modup.go
circuits/ckks/bootstrapping/fast_packing.go
circuits/ckks/bootstrapping/fast_bootstrap.go
```

Tests/benchmarks:

```text
circuits/ckks/bootstrapping/fast_bootstrap_bench_test.go
relevant Fast unit tests
```

Avoid broad unrelated refactors.

---

## Completion

Use a commit message reflecting the actual result. Preferred when production optimization is made:

```text
perf(ckks/fast): reduce Fast Bootstrap allocation overhead
```

If profiling proves no safe material optimization is worthwhile and only evidence/docs/benchmarks change:

```text
perf(ckks/fast): finalize Fast Bootstrap profiling
```

Push to:

```text
origin/fast-ckks
```

Final report must include:

1. exact profiling/benchmark commands;
2. Apple M4 + Go version/environment details available;
3. stable pre-optimization LogN13 Fast/Standard table;
4. stable pre-optimization LogN16 Fast/Standard table;
5. top CPU hotspots before optimization;
6. top alloc_space hotspots before optimization;
7. high-level stage attribution;
8. optimization(s) selected and why;
9. files changed;
10. whether compact dormant-row storage was implemented;
11. before/after Fast ns/op;
12. before/after Fast B/op;
13. before/after Fast allocs/op;
14. final Fast/Standard comparison;
15. repeated warm Fast benchmark result;
16. post-optimization CPU/alloc profile summary;
17. remaining dominant bottleneck;
18. correctness/oracle regression results;
19. confirmation no Standard/QP/Gadget/key-switch fallback;
20. any architecture doc update;
21. commit hash;
22. push result;
23. worktree/remote status.
