# Fast-CKKS Q-prefix repository guidance

> Instruction revision: `FAST-QPREFIX-ACTIVE-001`
>
> If a running agent claims this file requires `fast-ckks`, that agent is using stale or different instructions. On this revision, the required working branch is `fast-qprefix`.

## Active development line

The authoritative active development branch for the current Q-prefix work is:

`fast-qprefix`

The historical `fast-ckks` branch is a separate predecessor/comparison line and must not be used as the working branch for Q-prefix tasks unless an explicit future task says otherwise.

For architecture:
- `docs/FAST_QPREFIX_SPEC.md` is authoritative for Q-prefix branch architecture and policy.
- `docs/FAST_CKKS_SPEC.md` is inherited historical Fast-CKKS background. Its branch-specific `fast-ckks` wording and old implementation-status sections are not authoritative for `fast-qprefix` when they conflict with `FAST_QPREFIX_SPEC.md`, `CURRENT_TASK.md`, or current source.

## Start here

Before starting any task on this branch:

1. Safely synchronize the local `fast-qprefix` branch with `origin/fast-qprefix` before reading the task.
   - Check that the current branch is `fast-qprefix` and inspect `git status --short`.
   - If the worktree is clean, run `git fetch origin` and update with `git pull --ff-only origin fast-qprefix`.
   - If the worktree has uncommitted changes, the branch is not `fast-qprefix`, or the fast-forward update fails, do not reset, stash, overwrite, or discard anything automatically. Stop and report the condition instead.
   - After synchronization, re-read this `AGENTS.md` because the repository instructions themselves may have changed.
2. Read `CURRENT_TASK.md`.
3. Read the task specification referenced there.
4. Read `docs/FAST_QPREFIX_SPEC.md` for the authoritative Q-prefix architecture and invariants.
5. Read `docs/FAST_CKKS_SPEC.md` only for inherited Fast-CKKS background/invariants that do not conflict with the Q-prefix architecture.
6. Inspect the relevant current source before editing; repository evidence takes precedence over assumptions from names or prior knowledge.
7. Implement only the requested scope. Do not begin the next task implicitly.
8. Run the relevant targeted tests and benchmarks.
9. Review the diff and remove unrelated changes.
10. Commit and push the completed task to `origin/fast-qprefix` unless the task explicitly says otherwise. If push fails, keep the local commit intact and report the failure; do not rewrite history merely to retry.
11. Report the result, tests, benchmark changes, commit hash, and push result.

`CURRENT_TASK.md` is only the current work pointer. Task-specific requirements belong in the referenced Primary specification. Durable Q-prefix architecture belongs in `docs/FAST_QPREFIX_SPEC.md`. Do not duplicate those layers unnecessarily.

## Purpose

This `fast-qprefix` branch is the active intentionally insecure Fast-CKKS Q-prefix backend built on Lattigo. Keep application CKKS usage unchanged and keep Normal Lattigo behavior available. The current engineering priority is speed while preserving the bounded Q-prefix correctness contract.

The separate `fast-ckks` branch remains a historical/private-F comparison branch. Do not switch to it, merge it into this line, or rewrite it merely because older documentation mentions that branch.

## Persistent invariants

- `q_i` names the modulus at RNS index `i`; `r_i[k]` names coefficient `k`'s stored residue modulo `q_i`; `c0`, `c1`, `c2`, ... name ciphertext polynomial components.
- Q-prefix v2 maintains actual logical Q residues under the authoritative policy in `docs/FAST_QPREFIX_SPEC.md`; do not introduce private-F storage on this branch.
- Fast may maintain only the minimum/capped residue subset authorized by the Q-prefix architecture. Do not perform unnecessary full-RNS work, but follow Standard CKKS Level/Scale semantics when an operation requires level transitions, including transitions to Level 0.
- Never call Standard full-RNS code on stale or unmaintained Fast residue storage as a hidden fallback.
- The current zero-secret implementation is a mode, not a permanent scientific invariant. Preserve extension points for sampled secrets and future noise experiments.
- Preserve public APIs and CKKS parameter objects where practical; preserve structure while eliding expensive dormant or security-only computation.
- Do not add expensive noise-fidelity work during Stage A. Normal and Fast implementations must coexist for correctness and performance comparison.
- Keep changes CKKS-specific where possible. Do not modify BFV, BGV, TFHE, or unrelated infrastructure merely for consistency.

## Conflict resolution

If instructions appear inconsistent, use this precedence for the current Q-prefix line:

1. explicit current task specification;
2. `docs/FAST_QPREFIX_SPEC.md`;
3. current source and accepted test evidence;
4. inherited non-conflicting invariants from `docs/FAST_CKKS_SPEC.md`.

Do not silently invent a reconciliation. If a genuine semantic/architectural conflict remains after applying this precedence, stop and report it.

A stale mention of `fast-ckks` in historical documentation is not, by itself, a reason to leave `fast-qprefix`.

## Working style

Trace the exact source call path, identify the smallest modification surface, change one bounded behavior at a time, and prefer source-backed behavior over guesswork. Use the reusable fastdiag framework for supported performance attribution instead of rebuilding one-off instrumentation.
