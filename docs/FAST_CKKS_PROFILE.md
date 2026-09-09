# Fast Bootstrap profiling record

Task 010 profiling was run on an Apple M4, `darwin/arm64`, with Go
`go1.26.4`. The benchmark used the accepted Task 009-P profile:

```text
LogSlots=4, input Level=0, output Level=1
CosDiscrete, K=16, Mod1 degree=30, DoubleAngle=3, Mod1 inverse=0
```

The benchmark constructs parameters/evaluators/keys outside the timed child
benchmark and copies the input inside each timed iteration for both Fast and
Standard.

## Commands

Stable comparison before the optimization:

```bash
GOCACHE=/private/tmp/lattigo-go-cache go test ./circuits/ckks/bootstrapping \
  -run '^$' -bench 'BenchmarkFastBootstrapEndToEnd' -benchmem \
  -benchtime=3x -count=3
```

Fast-only CPU and allocation profiles used 200 steady-state iterations:

```bash
GOCACHE=/private/tmp/lattigo-go-cache go test ./circuits/ckks/bootstrapping \
  -run '^$' -bench 'BenchmarkFastBootstrapEndToEnd/LogN=13/Fast$' \
  -benchmem -benchtime=200x \
  -cpuprofile /private/tmp/fast-bootstrap-cpu-200.out \
  -memprofile /private/tmp/fast-bootstrap-mem-200.out
```

The Go 1.26 installation did not expose `go tool pprof`, so the generated
profiles were read with the cached `github.com/google/pprof/profile` parser.
No binary profile was added to the repository.

## Stable comparison before optimization

Values below are medians of the three `-count=3` runs. Standard is measured
with the same parameters and public API benchmark path.

| LogN | Path | ns/op | B/op | allocs/op |
|---:|---|---:|---:|---:|
| 13 | Fast | 10,479,653 | 11,312,194 | 3,437 |
| 13 | Standard | 182,514,722 | 37,380,170 | 21,131 |
| 16 | Fast | 89,978,778 | 98,239,077 | 3,772 |
| 16 | Standard | 1,946,114,708 | 295,556,410 | 23,584 |

## CPU profile before optimization

The 200-iteration LogN13 profile showed these dominant cumulative paths:

```text
Fast Bootstrap/bootstrapCore       ~1.52 s
Fast EvalMod                       ~0.94 s
Fast Rescale/rescaleN              ~0.92 s
Fast polynomial Evaluate            ~0.74 s
Fast DFT dft                        ~0.43 s
q0/q1 NTT                           ~0.39 s
runtime malloc/makeslice             ~0.85 s combined
```

The arithmetic entries are required Stage-A work. The runtime allocation
entries were the avoidable part investigated in this task.

## Allocation profile before optimization

The same profile showed full logical-level allocation as the material waste:

```text
ring.NewPoly / ring.Ring.NewPoly       ~1.84 GB / ~1.38 GB cumulative
rlwe.NewElement/NewCiphertext          ~1.07 GB cumulative
Fast bootstrapCore                      ~2.15 GB cumulative
Fast EvalMod                            ~1.28 GB cumulative
```

The source review confirmed that Fast-only constructors still created every
`q0...qLevel` backing array even though Fast reads/writes only q0/q1.

## Selected optimization

Added Fast-local compact ciphertext construction and resizing:

```text
logical Coeffs length = Level+1
Coeffs[0] and Coeffs[1] = N-sized maintained rows
Coeffs[2:] = nil dormant rows
```

Fast arithmetic, DFT, polynomial, Mod1, packing, ScaleDown, ModUp, and
rescale paths now use the compact helper where they own internal state.
Standard `ring.NewPoly`, `rlwe.NewElement`, `rlwe.NewCiphertext`, and all
Standard evaluator behavior remain unchanged. Public residual Level-0/1
outputs retain fully materialized active q0/q1 rows. Compact objects are not
sent to Standard/full-RNS APIs.

## Final Fast results

The final comparison used the same `-benchtime=3x -count=3` command. Values are
medians of the three runs.

| LogN | Path | ns/op | B/op | allocs/op |
|---:|---|---:|---:|---:|
| 13 | Fast | 9,692,250 | 2,405,538 | 3,331 |
| 13 | Standard | 177,232,375 | 37,439,744 | 21,131 |
| 16 | Fast | 88,147,083 | 19,078,389 | 3,655 |
| 16 | Standard | 1,962,251,361 | 294,684,373 | 23,581 |

Relative to the pre-optimization Fast medians, Fast B/op fell by about 79% at
LogN13 and 81% at LogN16. Fast time fell by about 7.5% and 2.0%, respectively.

## Warm repeated Fast result

```bash
GOCACHE=/private/tmp/lattigo-go-cache go test ./circuits/ckks/bootstrapping \
  -run '^$' -bench 'BenchmarkFastBootstrapEndToEnd/LogN=(13|16)/Fast$' \
  -benchmem -benchtime=10x -count=3
```

Median warm results:

| LogN | ns/op | B/op | allocs/op |
|---:|---:|---:|---:|
| 13 | 9,853,488 | 2,405,178 | 3,330 |
| 16 | 89,882,429 | 19,078,012 | 3,654 |

## Post-optimization profile and remaining bottleneck

The post-optimization 200-iteration LogN13 profile reported approximately
`9.93 ms/op`, `2.41 MB/op`, and `3330 allocs/op`. Full logical-level
`rlwe.NewCiphertext` allocation no longer appeared as the steady-state Fast
source of memory. The remaining per-call allocation is primarily the required
q0/q1 temporary state from `fast.NewCiphertext`/`compactPoly`; the remaining
CPU is dominated by EvalMod, polynomial multiplication/rescale, q0/q1 NTT, and
DFT mathematics. Reusing returned public ciphertext storage would violate
ownership, and concurrent scratch reuse is out of scope, so no further
speculative optimization was retained.

Approximate end-to-end stage attribution from the call stacks and existing
stage benchmarks is:

```text
EvalMod + polynomial/rescale       ~55-60% (overlapping cumulative estimate)
CoeffsToSlots + SlotsToCoeffs      ~20-25%
ModUp + Trace + packing/finalize   ~15-20%
```

These percentages are estimates, not additive instrumentation measurements.

For nearby stage evidence, the matching LogN13 microbenchmarks reported
`BenchmarkFastEvalModLogN13` at approximately `6.80 ms/op`, `1.36 MB/op`,
`2964 allocs/op`, and `BenchmarkFastPackingLogN13Count2` at approximately
`57.3 us/op`, `527 KB/op`, `30 allocs/op`. This is consistent with EvalMod
being the dominant isolated stage and packing being a small fraction of the
end-to-end path.

## Restrictions and correctness

The compact representation is Fast-local only. No new Standard evaluator,
full-Q/QP arithmetic, GadgetProduct, evaluation-key/key-switch fallback,
arbitrary-precision CRT, or parallel scratch reuse was introduced. Full-slot,
Level-1, BootstrapMany order/ownership, zero-secret ordinary decrypt/decode,
DFT planning parity, Standard-vs-Fast decode, and dormant-residue poison tests
remain required regressions.
