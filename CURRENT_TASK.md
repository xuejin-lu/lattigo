# Current Task

Task: QPREFIX-PERF-DIAG-007
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/QPREFIX-PERF-DIAG-007-POST-OPT002-RESCALE-REPROFILE.md`

Task class:
`D — Performance Diagnosis`

Purpose:
Freshly profile the current rows4 Q-prefix Rescale after OPT-002.

No Secondary production arithmetic changes.

Measure:
- diagnostics-off full rows4 Rescale;
- fresh CPU profile;
- refreshed fixed-width / residue / source-transform / restore-transform phase benchmarks;
- remaining division helper costs.

Do not:
- reuse DIAG-006 profile shares as current attribution;
- add per-coefficient timers;
- change Q-prefix policy/width;
- change transactionality;
- introduce reciprocal approximation, unsafe, assembly, or F/full-RNS fallback.

Return the required reprofile classification/metrics and `READY_FOR_WEB_REVIEW`.
