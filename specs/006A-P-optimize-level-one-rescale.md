# Task 006A-P — Optimize Fast Level-1 → Level-0 Rescale

## Goal

Remove the unnecessary CRT cost from the common Fast Bootstrap ScaleDown boundary:

```text
Level 1 -> Level 0
```

Task 006A correctness is already accepted. This task is performance-only and must preserve the existing higher-Level Fast Rescale architecture.

## Requirements

- In Fast `rescaleN`, specialize only when:

  ```text
  op0.Level() == 1
  nbRescales == 1
  ```

- For ordinary/non-Montgomery NTT representation, use the q0/q1-only Ring rounded-division primitive equivalent to:

  ```go
  ringQ.AtLevel(1).DivRoundByLastModulusNTT(...)
  ```

- This is permitted because Level 1 has valid maintained `r0` and `r1`, and the dropped divisor is exactly `q1`.

- Do not invoke Standard `ckks.Evaluator`, QP arithmetic, GadgetProduct, key switching, or full-Q execution.

- Preserve the existing fixed-width CRT path unchanged for:
  - Level >= 2; or
  - multi-level rescaling.

- Inspect source behavior before deciding Montgomery handling. If the Ring primitive is not valid for Montgomery-form input, leave Montgomery on the existing generic Fast path rather than adding conversion overhead.

- Support in-place and out-of-place operation.

- Preserve exact Standard rounded-rescale semantics:

  \[
  y=\operatorname{round}(x/q_1)
  \]

- Preserve metadata:

  ```text
  Level = 0
  Scale' = Scale / q1
  ```

## Tests

Add focused coverage for the specialized path:

- Level 1 -> 0
- ordinary NTT
- in-place
- out-of-place
- positive and negative coefficients
- rounding-boundary coefficients
- nonzero quotients
- exact `r0` comparison against Standard
- Scale and Level metadata

Keep all existing higher-Level and Montgomery Fast Rescale tests passing.

Also rerun:

```bash
go test ./schemes/ckks/fast
go test ./circuits/ckks/bootstrapping
go test ./circuits/ckks/dft
```

## Benchmark

Report Fast and Standard for:

1. direct Level 1 -> 0 Rescale;
2. Bootstrap ScaleDown at LogN=13;
3. Bootstrap ScaleDown at LogN=16.

Report:

```text
ns/op
B/op
allocs/op
```

Current ScaleDown baseline:

```text
LogN13
Fast      434,907 ns/op
Standard  197,614 ns/op

LogN16
Fast      3,668,765 ns/op
Standard  1,497,376 ns/op
```

The goal is to remove the current unnecessary ~2.4x penalty. A forced speedup over Standard is not required if both resolve to essentially the same two-limb Ring operation.

## Optional small cleanup

If trivial, compute `MulIntegerMaintained` scalar residues once per maintained limb rather than once per ciphertext component.

Do not broaden this into unrelated optimization.

## Out of scope

Do not implement:

- ModUp
- Trace
- sparse physical ciphertext storage
- QP paths
- full Bootstrap orchestration
- unrelated evaluator cleanup

## Preferred files

```text
schemes/ckks/fast/rescale.go
schemes/ckks/fast/rescale_test.go
schemes/ckks/fast/evaluator.go
circuits/ckks/bootstrapping/fast_scaledown_test.go
```

Do not modify Standard behavior.

## Completion

Commit with:

```text
perf(ckks/fast): specialize level-one Rescale
```

Push `origin/fast-ckks`.

Final report should include:

1. specialization condition;
2. Ring primitive used;
3. Montgomery behavior;
4. confirmation higher-Level CRT path is unchanged;
5. Level1->0 Fast/Standard benchmark;
6. ScaleDown LogN13 benchmark;
7. ScaleDown LogN16 benchmark;
8. B/op and allocs/op;
9. tests;
10. commit hash;
11. push result.
