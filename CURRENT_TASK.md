# Current Task

Task: QPREFIX-IMPL-009
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-IMPL-009-PERFORMANCE-RELEASE-GATE.md`

Production candidate:
`f9c7f21e65915bd3eafcd5b12590b570c22a7d6f`

Pre-F comparison point:
`40532b4dce5c7eeae2db5b0b6f21be64801ce923`

This is validation/release only.

Do not change production source.

Execute:
- detached baseline worktree;
- matched same-code benchmarks;
- identical temporary q0=55 P93 benchmark on both commits;
- current q0=56 P93 Fast + Standard comparison;
- full correctness/capacity/structural/fallback/public-output audit;
- final performance table and release verdict.

Keep q0=55 matched-baseline results separate from q0=56 production results.

Return exactly one release status:
- `QPREFIX_V2_RELEASE_PASS`
- `QPREFIX_V2_PERFORMANCE_REVIEW`
- `QPREFIX_V2_RELEASE_FAIL`

Then report `READY_FOR_WEB_REVIEW`.
