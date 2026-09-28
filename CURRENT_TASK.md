# Current Task

Task: QPREFIX-PERF-DIAG-003
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-PERF-DIAG-003-GENERATED-POWER-RESCALE-CAUSALITY.md`

Historical q0=55 baseline:
`40532b4dce5c7eeae2db5b0b6f21be64801ce923`

Q-prefix-v2 production candidate:
`f9c7f21e65915bd3eafcd5b12590b570c22a7d6f`

This is diagnosis only.

Use immutable detached comparison points for the matched q0=55 EvalMod-real generated-power phase.

Measure:
- actual generated-power runtime DAG and operation counts;
- per-power timing for {2,3,4,6,8,16};
- matched baseline/candidate Rescale cost;
- candidate Rescale preflight versus materialization;
- non-overlapping generated-power delta accounting.

Do not modify production code.
Do not remove or bypass preflight.
Do not alter transactional failure-before-mutation semantics.
Do not alter parameters, schedules, generated powers, or `QPrefixWidth(Level)`.
Do not introduce F/full-RNS fallback.

Return the required diagnostic classification and `READY_FOR_WEB_REVIEW`.
