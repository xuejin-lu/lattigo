# Current Task

Task: QPREFIX-PERF-DIAG-004
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/QPREFIX-PERF-DIAG-004-POST-OPT-HOTSPOT-REFRESH.md`

Task class:
`D — Performance Diagnosis`

Purpose:
Measure the post-optimization hotspot distribution after transactional Rescale staging.

No Secondary production modification is allowed.

Use the reusable fastdiag framework only:
- current stage trace;
- current stage,power,rescale trace;
- compare pre-opt `74cb73dcff6c552cda0671ed7faea897b448fbbd` against post-opt `d50ff4db757d4a2b9922937a4e7f316fd3f286b9`.

Do not:
- optimize Rescale or another kernel;
- add one-off overlays;
- alter Q-prefix policy, parameters, schedules, or F/full-RNS architecture.

Return the required diagnostic result and `READY_FOR_WEB_REVIEW`.
