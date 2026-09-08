# Fast-CKKS repository guidance

## Start here

Before starting any Fast-CKKS task:

1. Safely synchronize the local `fast-ckks` branch with `origin/fast-ckks` before reading the task.
   - Check that the current branch is `fast-ckks` and inspect `git status --short`.
   - If the worktree is clean, run `git fetch origin` and update with `git pull --ff-only origin fast-ckks`.
   - If the worktree has uncommitted changes, the branch is not `fast-ckks`, or the fast-forward update fails, do not reset, stash, overwrite, or discard anything automatically. Stop and report the condition instead.
   - After synchronization, re-read this `AGENTS.md` because the repository instructions themselves may have changed.
2. Read `CURRENT_TASK.md`.
3. Read the task specification referenced there.
4. Read `docs/FAST_CKKS_SPEC.md` for durable architecture and invariants.
5. Inspect the relevant current source before editing; repository evidence takes precedence over assumptions from names or prior knowledge.
6. Implement only the requested scope. Do not begin the next task implicitly.
7. Run the relevant targeted tests and benchmarks.
8. Review the diff and remove unrelated changes.
9. Commit and push the completed task to `origin/fast-ckks` unless the task explicitly says otherwise. If push fails, keep the local commit intact and report the failure; do not rewrite history merely to retry.
10. Report the result, tests, benchmark changes, commit hash, and push result.

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
