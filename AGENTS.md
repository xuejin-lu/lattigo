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
- `docs/FAST_CKKS_SPEC.md` is historical Fast-CKKS background, not part of the normal startup reading chain and not an authority for current Q-prefix behavior. Consult it only when the active Primary task explicitly requires historical context.

## Start here

Before starting any task on this branch:

1. Safely synchronize the local `fast-qprefix` branch with `origin/fast-qprefix` before reading the task.
   - Check that the current branch is `fast-qprefix` and inspect `git status --short`.
   - If the worktree is clean, run `git fetch origin` and update with `git pull --ff-only origin fast-qprefix`.
   - If the worktree has uncommitted changes, the branch is not `fast-qprefix`, or the fast-forward update fails, do not reset, stash, overwrite, or discard anything automatically. Stop and report the condition instead.
   - After synchronization, re-read this `AGENTS.md` because the repository instructions themselves may have changed.
2. Safely synchronize the Primary `xuejin-lu/heart-lattigo-bootstrap` `main` branch according to its `AGENTS.md`; read its freshly synchronized `CURRENT_TASK.md` and the active specification referenced there. A blocked, completed, or missing Primary task does not authorize Secondary execution.
3. Read `docs/FAST_QPREFIX_SPEC.md` for authoritative Q-prefix architecture and invariants.
4. Inspect the relevant current source before editing; repository evidence takes precedence over assumptions from names or prior knowledge.
5. Consult historical `docs/FAST_CKKS_SPEC.md` only if the active Primary task explicitly needs historical design context. Do not infer active instructions from this repository's `CURRENT_TASK.md`, which is only a pointer to the Primary task authority.
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

Normative rules and factual evidence are different and must not share one precedence list.

For the current Q-prefix line, use this **normative precedence**:

1. durable repository invariants and prohibitions in this `AGENTS.md`, together with the authoritative Q-prefix mathematical/architectural contract in `docs/FAST_QPREFIX_SPEC.md`;
2. any explicit current experiment/architecture contract that has already been accepted without changing the durable rules above;
3. the current task specification;
4. `CURRENT_TASK.md`, which is only a pointer to the active task and must not introduce new architecture or semantics;
5. inherited historical material from `docs/FAST_CKKS_SPEC.md` only where it is explicitly non-conflicting.

A task specification must not silently override a durable invariant or the authoritative Q-prefix mathematics. If a task requires behavior that conflicts with a higher-level rule, stop and report `NEEDS_WEB_REVIEW`; do not choose the newer task merely because it is more recent.

Current source code, tests, benchmark outputs, result artifacts, and Git history are **evidence about what the implementation currently does**, not normative authority about what it is allowed to do. If source or accepted tests conflict with a durable rule or task contract, report the mismatch as an implementation/specification defect; do not use implementation reality to redefine the rule.

If new research genuinely requires changing a durable invariant or Q-prefix architectural rule, treat that as an explicit constitution/architecture amendment:
1. state the rule being changed and the evidence/rationale;
2. obtain Web scientific/architectural review;
3. update the authoritative durable document first;
4. only then create an implementation task under the amended rule.

Do not create ad-hoc per-task exceptions to bypass this amendment process.

A stale mention of `fast-ckks` in historical documentation is not, by itself, a reason to leave `fast-qprefix`.

## Working style

Trace the exact source call path, identify the smallest modification surface, change one bounded behavior at a time, and prefer source-backed behavior over guesswork. Use the reusable fastdiag framework for supported performance attribution instead of rebuilding one-off instrumentation.
