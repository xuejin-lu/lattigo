# Current Task

Task: QPREFIX-PERF-DIAG-005
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/QPREFIX-PERF-DIAG-005-LOW-OVERHEAD-POWER-ATTRIBUTION.md`

Task class:
`D — Performance Diagnosis`

Purpose:
Validate the post-opt generated-power bottleneck with low instrumentation perturbation.

Use only:
- `--trace power`
- `--trace stage,power`

Do NOT enable `rescale` scope for the authoritative category ranking.

No Secondary production changes.
No ad-hoc overlays.
Do not optimize any kernel.

Return the required diagnostic classification/metrics and `READY_FOR_WEB_REVIEW`.
