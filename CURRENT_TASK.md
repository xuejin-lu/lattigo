# Current Task

Task: FAST-STANDARD-NUMERICAL-FIX-001
Status: READY_FOR_CODEX

Working branch:
`fast-qprefix`

Required marker:
`FAST-QPREFIX-ACTIVE-001`

This task **supersedes and suspends**:
`QPREFIX-PERF-MEASURE-LOGN16-001`

Do not resume the LogN16 measurement-only task unless a later explicit task restores it.

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap/specs/FAST-STANDARD-NUMERICAL-FIX-001-REMOVE-NORMALIZED-EVALMOD.md`

Production parent before this repair:
`97c1c6174e0d7ef781de6d8dadce5c2869a53496`

Goal:
Remove the historical Fast-only normalized LogN13 EvalMod production path and make Fast Q-prefix EvalMod follow genuine Standard Mod1 mathematics and Scale/Level progression while retaining Fast Q-prefix storage/arithmetic.

Authorized scope:
- Secondary production code changes required by FIX-001 are explicitly authorized.
- Modify only the bounded EvalMod/Fast-Q-prefix surfaces required by the Primary spec.
- Preserve Standard production arithmetic.
- Preserve the current zero-a / error-retaining key semantics.
- Do not reintroduce Standard full-RNS/QP/key-switch fallback.
- Use Q0123 capacity guards as assertions; if the Standard-equivalent schedule actually fails capacity, stop and report the exact first deficit instead of inventing another workaround.

The prior DIAG-003 semantic-alignment path is superseded.
Do not restore or regenerate its abandoned result files.
Do not tune the historical normalized LogN13 exponents.

Required Secondary workflow:
1. synchronize `fast-qprefix`;
2. re-read `AGENTS.md`, this file, the authoritative Primary FIX-001 spec, and `docs/FAST_QPREFIX_SPEC.md`;
3. implement FIX-001;
4. run the required Secondary tests;
5. review the diff for unrelated changes;
6. commit and push the authorized Secondary production repair to `origin/fast-qprefix`;
7. use the resulting Secondary commit for the Primary canonical numerical rerun and FIX-001 result artifact.

Required completion is governed by the Primary FIX-001 spec:
`FAST_STANDARD_EVALMOD_FIX_READY`
then
`READY_FOR_WEB_REVIEW`.
