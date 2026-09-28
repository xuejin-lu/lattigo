# Current Task

Task: QPREFIX-PERF-DIAG-002
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-PERF-DIAG-002-MATCHED-Q55-STAGE-CLOSURE.md`

Historical baseline:
`40532b4dce5c7eeae2db5b0b6f21be64801ce923`

Q-prefix-v2 production candidate:
`f9c7f21e65915bd3eafcd5b12590b570c22a7d6f`

This is diagnosis only.

Important:
- use immutable detached comparison points for the matched q0=55 replay;
- do not move or rewrite the authoritative `fast-qprefix` worktree for measurement;
- measure the exact same q0=55 P93 workload on baseline and candidate;
- capture non-overlapping in-context Bootstrap stage timings;
- identify the single largest positive stage delta and drill down only that stage;
- no production optimization or architecture change is authorized.

Do not modify production code.
Do not alter parameters, schedules, or `QPrefixWidth(Level)`.
Do not introduce F/full-RNS fallback.

Return the required diagnostic status and `READY_FOR_WEB_REVIEW`.
