# Fast CKKS — Phase 1 Constitution

**Project:** Lattigo Fast CKKS  
**Branch:** `fast-ckks`  
**Current checkpoint:** `c190fb624d2f3e78dda3b32a71ddd067dd67d0d7`  
**Status:** Phase 1C Step 3C committed; later stages not started.

---

# 1. Purpose

This document is the architectural constitution for the Fast CKKS implementation.

It records the mathematical assumptions, execution invariants, design decisions, scope boundaries, and implementation sequence that have been established for Phase 1.

Future implementation and review must follow this document.

A future implementation agent must **not reinterpret or silently replace these decisions with Standard Lattigo semantics** merely because the Standard implementation follows a different design.

If a future stage appears to contradict this constitution, the correct action is:

1. stop implementation;
2. identify the contradiction;
3. report the architectural conflict;
4. request an explicit architecture decision.

Do not silently change the Fast semantics.

---

# 2. Core Design Goal

The Fast CKKS implementation intentionally sacrifices the security properties of the standard RLWE construction in exchange for computational simplification.

This is intentional.

The goal is **not** to preserve Standard Lattigo's cryptographic security model.

The goal is to construct a computationally simplified CKKS execution path suitable for the intended research/demo environment.

The central simplification is:

```text
Secret key:

    s = 0
```

Consequently, the key-dependent deterministic term disappears.

This is the fundamental architectural assumption underlying the Fast path.

---

# 3. Fundamental Fast-Key Semantics

Fast keys use the structure:

```text
(error, 0)
```

instead of the Standard RLWE-style:

```text
(error + deterministic term, a)
```

The second component is explicitly zero.

The important consequence is that multiplication by the secret key becomes trivial:

```text
a * s = a * 0 = 0
```

and therefore the deterministic secret-key-dependent component disappears.

This is intentional and is the reason the Fast implementation can remove substantial amounts of Standard RLWE computation.

---

# 4. The `s = 0` Decision

The following is an explicit architecture decision:

```text
FAST SECRET KEY:

    s = 0
```

This is **not** an optimization heuristic.

It is part of the mathematical definition of the Fast execution model.

Do not reintroduce a non-zero secret key into Fast execution merely to preserve Standard RLWE semantics.

---

# 5. Degree-2 Ciphertexts and `c₂`

A degree-2 ciphertext has the conceptual form:

```text
c0 + c1*s + c2*s²
```

Under the Fast assumption:

```text
s = 0
```

therefore:

```text
c1*s = 0
c2*s² = 0
```

and the expression reduces to:

```text
c0
```

More generally, the higher-degree secret-key terms do not contribute to the Fast decryption semantics.

Therefore the Fast multiplication path is explicitly allowed to produce:

```text
degree 2
```

temporarily and then immediately discard the `c2` component.

The Fast operation is:

```text
degree-1 × degree-1
        ↓
degree-2 temporary result
        ↓
discard c2
        ↓
degree-1 result
```

This is **not Standard relinearization**.

It must not call:

```text
GadgetProduct
```

and it must not generate a Standard relinearization correction term.

---

# 6. Fast Relinearization Semantics

There is no Standard-style Fast relinearization.

The Fast path does not attempt to transform:

```text
c2*s² → c2*s
```

because:

```text
s = 0
```

makes the term irrelevant to the Fast semantics.

Therefore:

```text
Fast multiplication
    = polynomial multiplication of c0/c1
      followed by direct degree truncation
```

not:

```text
Standard multiplication
    + Standard relinearization
```

`FastTruncateDegree2To1` is therefore a valid Fast operation.

---

# 7. Authoritative RNS Limbs

The most important execution invariant is:

```text
q0, q1 = authoritative
q2 ... qL = dormant / stale / non-authoritative
```

Only `q0` and `q1` are required to represent the current Fast computation.

The Fast execution path therefore does **not** maintain mathematically synchronized values in:

```text
q2 ... qL
```

This is intentional.

---

# 8. Why Only q0/q1 Matter

The Fast architecture intentionally uses the first two RNS components as the authoritative representation.

The current Fast design relies on the fact that two sufficiently large RNS components can represent/reconstruct the required coefficient information for the intended Fast arithmetic.

The remaining RNS components are therefore not part of the authoritative computational state.

They are storage that may contain:

```text
old values
stale values
constructor-initialized values
```

and these values are not considered part of the current mathematical Fast state.

Do not infer correctness of q2...qL from their stored values.

---

# 9. Dormant-Limb Rule

For Fast production execution:

```text
q2 ... qL:

    MUST NOT be treated as authoritative.
```

Unless a future architecture decision explicitly changes this rule, Fast operations should:

- not read q2...qL;
- not multiply q2...qL;
- not perform NTT on q2...qL;
- not perform INTT on q2...qL;
- not synchronize q2...qL after every operation;
- not reconstruct q2...qL after every multiplication;
- not perform full-RNS arithmetic merely to keep dormant limbs current.

If a Standard operation requires those limbs, that is an explicit Fast/Standard boundary problem and must be handled by a deliberately designed boundary operation.

It must not cause every Fast operation to maintain full-RNS state.

---

# 10. No Per-Operation CRT

CRT reconstruction is **not part of the Fast hot path**.

The old reference implementation may contain:

```text
ReconstructQ0Q1
CheckCoefficientBounds
Redistribute
```

These are reference/materialization mechanisms.

They must not automatically execute after every Fast multiplication.

The Fast production path is:

```text
q0/q1
    ↓
Fast NTT
    ↓
q0/q1 arithmetic
    ↓
Fast INTT
    ↓
q0/q1 authoritative result
```

No automatic:

```text
CRT
```

and no automatic:

```text
redistribution
```

after each operation.

---

# 11. No Per-Operation Redistribution

The following design is explicitly rejected:

```text
Fast multiplication
    ↓
CRT
    ↓
reconstruct x
    ↓
redistribute x into q0...qL
```

after every operation.

This would destroy the intended performance advantage.

The fact that q2...qL may be stale is intentional.

They do not need to be refreshed merely because q0/q1 changed.

---

# 12. NTT / INTT Domain Rule

Fast execution uses partial NTT/INTT.

The domain invariant is:

```text
q0/q1:
    coefficient domain ↔ NTT domain

q2...qL:
    remain in coefficient-domain storage
    and are dormant/non-authoritative
```

Therefore a valid Fast state can conceptually be:

```text
q0/q1:
    NTT domain

q2...qL:
    stale coefficient-domain storage
```

or:

```text
q0/q1:
    coefficient domain

q2...qL:
    stale coefficient-domain storage
```

The generic `ring.Poly` abstraction must not be modified merely to represent this Fast-specific mixed state.

Fast callers are responsible for maintaining the invariant.

---

# 13. Partial NTT

The Fast implementation must not call:

```text
ring.Ring.NTT
```

when doing Fast partial transforms because Standard `ring.Ring.NTT` processes all active RNS limbs.

Instead:

```text
ringQ.SubRings[0].NTT
ringQ.SubRings[1].NTT
```

are used.

Likewise for INTT:

```text
ringQ.SubRings[0].INTT
ringQ.SubRings[1].INTT
```

The current primitive is:

```text
FastPartialNTT
FastPartialINTT
```

These operate only on q0/q1.

---

# 14. Fast Multiplication

The authoritative Fast multiplication primitive must only operate on q0/q1.

Conceptually:

```text
A = (a0, a1)
B = (b0, b1)
```

Then:

```text
c0 = a0 * b0

c1 = a0 * b1 + a1 * b0

c2 = a1 * b1
```

All multiplication arithmetic is performed only on q0/q1.

The temporary degree-2 result is then truncated:

```text
(c0, c1, c2)
        ↓
(c0, c1)
```

No Standard relinearization is performed.

---

# 15. Fast Addition/Subtraction

Fast Add/Sub operates only on:

```text
q0/q1
```

and leaves:

```text
q2...qL
```

dormant.

The Standard full-RNS:

```text
ring.Ring.Add
ring.Ring.Sub
```

must not be used by Fast Add/Sub when the goal is q0/q1-authoritative execution.

Current primitives:

```text
FastAdd
FastSub
```

---

# 16. Fast Evaluator Boundary

Fast execution must be explicitly selected.

Current architecture:

```text
fast.Evaluator
```

rather than modifying:

```text
ckks.Evaluator
```

The Standard evaluator remains Standard.

Do not introduce:

```text
global FastMode
```

Do not detect Fast mode using coefficient values.

Do not infer Fast mode from:

```text
q0 == 0
q1 == 0
```

or any other numerical heuristic.

Do not add a generic:

```text
IsFast
```

flag to `rlwe.Ciphertext` unless a future explicit architecture decision approves such a change.

---

# 17. Standard/Fast Isolation

The intended architecture is:

```text
Standard:

ckks.Evaluator
    ↓
Standard RLWE operations
    ↓
full-RNS execution
```

versus:

```text
Fast:

fast.Evaluator
    ↓
q0/q1-authoritative primitives
    ↓
Fast execution semantics
```

The Standard path must remain functional and regression-tested.

Fast behavior must not silently activate for Standard ciphertexts.

---

# 18. Evaluation Keys

Fast evaluation keys are intentionally generated using the Fast key semantics.

Their relevant structure is:

```text
(error, 0)
```

Fast keys are not interchangeable with Standard evaluation keys.

Standard GadgetProduct must not silently consume Fast keys.

Fast key misuse should be rejected explicitly at established boundaries.

Fast key serialization/compression/expansion behavior is intentionally restricted where necessary.

---

# 19. Security Model

Security is **not** the optimization target of Phase 1.

The Fast architecture intentionally removes secret-key-dependent terms.

This means:

```text
security-preserving RLWE semantics
```

are not the objective of the Fast path.

Do not reject an optimization merely because it would be inappropriate for a secure RLWE implementation, provided that:

1. it is mathematically consistent with the Fast model;
2. it preserves the q0/q1 authoritative invariant;
3. it does not silently modify Standard execution;
4. it is within the approved phase scope.

---

# 20. Phase 1A — Fast RNS Reconstruction

Status:

```text
COMMITTED
```

Commit:

```text
3260b562 feat(ckks): implement Phase 1A Fast RNS reconstruction
```

Phase 1A established the reference reconstruction mechanisms:

```text
ReconstructQ0Q1
CheckCoefficientBounds
Redistribute
```

These remain useful as:

- correctness reference;
- debugging oracle;
- materialization boundary;
- test oracle.

They are not the intended Fast hot-path implementation.

---

# 21. Phase 1B — Fast Key Generation

Status:

```text
COMMITTED
```

Commit:

```text
5da28cb4 feat(ckks/fast): implement Fast CKKS key generation
```

Phase 1B established:

- Fast key layout;
- Fast PublicKey;
- Fast EvaluationKey;
- Fast Galois keys;
- Fast relinearization-key material;
- Bootstrap key-generation entry points;
- Fast-key misuse protection.

The key conceptual rule is:

```text
Fast key = (error, 0)
```

and the Fast secret-key semantics are:

```text
s = 0
```

---

# 22. Phase 1C Step 3A — Fast Add/Sub

Status:

```text
COMPLETED
```

Implemented:

```text
FastAdd
FastSub
```

Characteristics:

- q0/q1 only;
- no full-RNS arithmetic;
- q2...qL remain dormant;
- Standard evaluator unchanged.

---

# 23. Phase 1C Step 3B — Fast Evaluator / Multiplication Boundary

Status:

```text
COMMITTED
```

Commit:

```text
0a36f454 feat(ckks/fast): implement q0/q1-authoritative evaluator operations
```

Established:

```text
fast.Evaluator
```

and the authoritative multiplication boundary.

The Fast evaluator explicitly selects Fast execution.

---

# 24. Phase 1C Step 3C — Fast Degree Truncation

Status:

```text
COMMITTED
```

Commit:

```text
c190fb62 feat(ckks/fast): implement q0/q1-authoritative multiplication
```

Established:

```text
FastTruncateDegree2To1
```

Fast multiplication semantics are therefore:

```text
degree-1 × degree-1
    ↓
degree-2
    ↓
discard c2
    ↓
degree-1
```

No Standard relinearization.

---

# 25. Current Git Checkpoint

Current HEAD:

```text
c190fb62
```

Recent history:

```text
c190fb62 feat(ckks/fast): implement q0/q1-authoritative multiplication
0a36f454 feat(ckks/fast): implement q0/q1-authoritative evaluator operations
5da28cb4 feat(ckks/fast): implement Fast CKKS key generation
```

Earlier Phase 1A checkpoint:

```text
3260b562 feat(ckks): implement Phase 1A Fast RNS reconstruction
```

Current working tree at checkpoint:

```text
clean
```

Repository-wide tests:

```text
go test ./...
    PASS
```

---

# 26. Completed Implementation Inventory

The following Fast primitives currently exist:

```text
FastPartialNTT
FastPartialINTT

FastAdd
FastSub

FastMulQ01
FastMulQ01Authoritative

FastTruncateDegree2To1

fast.Evaluator
```

Reference mechanisms:

```text
ReconstructQ0Q1
CheckCoefficientBounds
Redistribute
```

---

# 27. Performance Lessons Already Established

An early benchmark of the reference Fast multiplication path showed that the dominant costs were:

```text
Redistribution
    >
CRT reconstruction
    >
bound checking / allocation
```

rather than:

```text
Partial NTT
q0/q1 multiplication
```

This is why the production authoritative path must avoid:

```text
CRT
redistribution
per-coefficient big.Int
```

during ordinary Fast execution.

Do not reintroduce these operations merely for convenience.

---

# 28. Fixed-Width CRT Investigation

A fixed-width CRT design was investigated.

The conclusion was:

- q0/q1 commonly use approximately 50–60-bit primes;
- q0*q1 is therefore commonly around 100–120 bits;
- a 128-bit equivalent representation such as `[2]uint64` is sufficient for current common parameters;
- fixed-width arithmetic can theoretically remove most per-coefficient allocations;
- CRT + redistribution can theoretically be fused.

However:

**this is not the current authoritative hot path.**

The reference big.Int implementation should remain available as an oracle.

Do not replace the reference implementation without an explicit task.

---

# 29. Important Architectural Consequence

The most important performance principle is:

> Do not repeatedly materialize information that Fast execution intentionally treats as non-authoritative.

In particular:

```text
Fast operation
    ≠
Fast operation + full-RNS synchronization
```

The latter defeats the purpose of the architecture.

---

# 30. Bootstrap Direction

The eventual goal includes Fast Bootstrap.

However, Standard Bootstrap cannot simply consume Fast ciphertexts because Standard Bootstrap assumes full-RNS semantics and Standard execution primitives.

The intended future direction is to replace relevant execution layers with Fast-specific primitives while preserving the high-level Bootstrap algorithm where mathematically compatible.

Potential future components include:

```text
Fast ModUp
Fast ModDown
Fast KeySwitch
Fast Rotation
Fast Rescale
Fast DFT
Fast EvalMod
Fast ring-switch operations
Fast Bootstrap orchestration
```

These must be implemented only after their mathematical contracts are explicitly established.

Do not prematurely modify the Standard Bootstrap evaluator.

---

# 31. Future Phase Sequence

The exact future sequence is subject to architecture review, but the current high-level dependency graph is:

```text
Phase 1A
    Fast RNS reconstruction/reference
        ↓
Phase 1B
    Fast key generation
        ↓
Phase 1C Step 3A
    Fast Add/Sub
        ↓
Phase 1C Step 3B
    Fast evaluator + authoritative multiplication
        ↓
Phase 1C Step 3C
    Fast degree truncation
        ↓
Future:
    Fast Relinearization semantics
    Fast KeySwitch
    Fast Rotation
    Fast Rescale
    Fast ModUp/ModDown
    Fast DFT
    Fast EvalMod
    Fast Bootstrap
```

The precise ordering of the future stages must be determined by dependency and mathematical necessity.

---

# 32. Rules for Future Codex Tasks

Every future Codex task should follow this sequence:

### Step 1 — Read this constitution

Do not start implementation before understanding the Fast invariants.

### Step 2 — Inspect the Standard source

Use Standard Lattigo as the reference for:

- call graph;
- metadata semantics;
- edge cases;
- level handling;
- alias handling;
- domain handling;
- mathematical reference.

But do **not** automatically copy Standard full-RNS computation into Fast.

### Step 3 — Identify the Fast mathematical contract

Explicitly state:

```text
What is authoritative?
What is dormant?
What must be computed?
What must NOT be computed?
What metadata must survive?
```

### Step 4 — Define the smallest implementation boundary

Avoid broad refactoring.

Prefer:

```text
Fast-specific primitive
```

or:

```text
Fast-specific evaluator entry point
```

over modifying generic infrastructure.

### Step 5 — Implement only the requested stage

Do not implement future stages opportunistically.

For example, while implementing Fast Rotation, do not silently implement:

```text
Fast Bootstrap
Fast Rescale
Fast EvalMod
```

### Step 6 — Test against Standard where appropriate

Correctness tests should compare the relevant mathematical result against the appropriate reference.

However, tests must respect the Fast semantics.

Do not require dormant q2...qL to match Standard values unless the explicit contract requires materialization.

### Step 7 — Audit forbidden work

Every Fast hot-path implementation should explicitly verify that it does not accidentally execute:

```text
full-RNS NTT
full-RNS INTT
full-RNS multiplication
per-operation CRT
per-operation redistribution
unnecessary q2...qL synchronization
```

### Step 8 — Benchmark

Measure:

```text
Fast
vs
Standard
```

but ensure the benchmark comparison is genuinely apples-to-apples.

Do not claim a speedup if one side is already in NTT domain and the other is not.

### Step 9 — Review before commit

No commit should be created until:

- correctness review passes;
- scope review passes;
- Standard regression passes;
- working tree scope is understood.

---

# 33. Explicitly Forbidden Shortcuts

Do not introduce:

```text
global FastMode
```

Do not use:

```text
coefficient-value heuristic
```

Do not use:

```text
q0/q1 == 0
```

as Fast-mode detection.

Do not silently modify:

```text
ring.Ring
```

to make Standard operations behave like Fast operations.

Do not modify:

```text
ring.Poly
```

to force Fast mixed-domain semantics globally.

Do not update dormant limbs simply because they exist in memory.

Do not insert CRT/redistribution after every Fast operation.

Do not call Standard full-RNS primitives from a Fast hot path merely because they are convenient.

Do not silently restore Standard relinearization semantics.

---

# 34. Standard Compatibility Rule

Standard Lattigo behavior must remain unchanged unless a future task explicitly states otherwise.

Every Fast optimization must be isolated from Standard execution.

After significant changes, run:

```text
go test ./...
```

and verify Standard packages remain healthy.

---

# 35. Architecture Escalation Rule

If an implementation agent encounters a situation where:

```text
Fast invariant
```

and:

```text
Standard API contract
```

cannot both be satisfied, the agent must **not invent a compromise**.

Instead report:

```text
ARCHITECTURE BLOCKER
```

and explain:

1. the Standard contract;
2. the Fast contract;
3. the conflict;
4. the minimum possible architecture change.

Only after an explicit architecture decision should implementation continue.

---

# 36. Current Non-Goals

The following are intentionally not part of the current completed implementation:

```text
Fast KeySwitch
Fast Rotation
Fast Rescale
Fast ModUp
Fast ModDown
Fast DFT
Fast EvalMod
Fast Bootstrap execution
```

They must be treated as future stages.

---

# 37. Final Constitution Summary

The Fast CKKS architecture is based on the following principles:

```text
1. Fast secret semantics:
       s = 0

2. Fast evaluation keys:
       (error, 0)

3. Fast ciphertext authority:
       q0/q1 authoritative

4. Dormant limbs:
       q2...qL stale/non-authoritative

5. Fast NTT:
       q0/q1 only

6. Fast INTT:
       q0/q1 only

7. Fast arithmetic:
       q0/q1 only

8. No routine full-RNS synchronization.

9. No per-operation CRT.

10. No per-operation redistribution.

11. Degree-2 multiplication:
       compute temporarily,
       discard c2,
       return degree 1.

12. No Standard GadgetProduct for Fast relinearization.

13. Fast execution is explicitly selected through
       fast.Evaluator.

14. No global FastMode.

15. No coefficient heuristics.

16. Standard Lattigo execution remains isolated.

17. When Fast and Standard contracts conflict,
       stop and request architecture review.

18. Optimize the authoritative state,
       not the dormant state.
```

---

# 38. One-Sentence Rule

The entire architecture can be summarized as:

> **Fast CKKS only computes and maintains what the Fast mathematical model actually considers authoritative: `s = 0`, q0/q1 are authoritative, q2...qL are dormant, and full-RNS materialization is an explicit boundary operation—not something performed after every Fast operation.**