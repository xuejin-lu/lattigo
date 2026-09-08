# Fast-CKKS repository guidance

## Start here

Before starting any Fast-CKKS task:

1. Read `CURRENT_TASK.md`.
2. Read the task specification referenced there.
3. Read `docs/FAST_CKKS_SPEC.md` for durable architecture and invariants.
4. Inspect the relevant current source before editing; repository evidence takes precedence over assumptions from names or prior knowledge.
5. Implement only the requested scope. Do not begin the next task implicitly.
6. Run the relevant targeted tests and benchmarks.
7. Review the diff and remove unrelated changes.
8. Report the result, tests, benchmark changes, commit hash, and push result.

`CURRENT_TASK.md` is only the current work pointer. Task-specific requirements belong in `specs/`. Durable Fast-CKKS architecture belongs in `docs/FAST_CKKS_SPEC.md`. Do not duplicate those layers unnecessarily.

## Purpose

This `fast-ckks` branch is an intentionally insecure Fast-CKKS backend built on Lattigo. Keep application CKKS usage unchanged and keep Normal Lattigo behavior available. The current engineering priority is speed: Fast execution must become substantially faster than Normal Lattigo for large workloads.

## Persistent invariants

- `q_i` names the modulus at RNS index `i`; `r_i[k]` names coefficient `k`'s stored residue modulo `q_i`; `c0`, `c1`, `c2`, ... name ciphertext polynomial components.
- Fast may maintain only the minimum residue subset needed for the current operation. Do not perform unnecessary full-RNS work, but follow Standard CKKS Level/Scale semantics when an operation requires level transitions, including transitions to Level 0.
- Never call Standard full-RNS code on stale or unmaintained Fast residue storage as a hidden fallback.
- The current zero-secret implementation is a mode, not a permanent scientific invariant. Preserve extension points for sampled secrets and future noise experiments.
- Preserve public APIs and CKKS parameter objects where practical; preserve structure while eliding expensive dormant or security-only computation.
- Do not add expensive noise-fidelity work during Stage A. Normal and Fast implementations must coexist for correctness and performance comparison.
- Keep changes CKKS-specific where possible. Do not modify BFV, BGV, TFHE, or unrelated infrastructure merely for consistency.

For detailed architecture, current implementation status, hot-path rules, validation, and future experimental knobs, use `docs/FAST_CKKS_SPEC.md`.

## Working style

Trace the exact source call path, identify the smallest modification surface, change one bounded behavior at a time, and prefer source-backed behavior over guesswork. If a task specification conflicts with the durable architecture specification or current source evidence, do not silently invent a reconciliation; report the conflict explicitly.
