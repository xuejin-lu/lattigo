# Current Task

Task: QPREFIX-PERF-DIAG-006
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/QPREFIX-PERF-DIAG-006-RESCALE-KERNEL-ATTRIBUTION.md`

Task class:
`D — Performance Diagnosis`

Purpose:
Attribute current rows4 Q-prefix Rescale cost with diagnostics-off focused benchmarks and CPU profiling.

Do not:
- modify production arithmetic;
- enable deep per-coefficient Rescale tracing for authoritative attribution;
- change Q-prefix policy/width;
- change P93/generated-power schedules;
- introduce F/full-RNS fallback.

Test/benchmark-only reusable diagnostics code is allowed if it does not affect production behavior.

Return the required attribution result and `READY_FOR_WEB_REVIEW`.
