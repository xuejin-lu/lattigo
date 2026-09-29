# Current Task

Task: FAST-STANDARD-NUMERICAL-001
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/FAST-STANDARD-NUMERICAL-001-P93-REFERENCE.md`

Task class:
`C — Numerical Correctness Validation`

Purpose:
Measure the actual P93 q0=55 numerical gap between current Fast Bootstrap and Standard Bootstrap on the canonical 4096-slot workload.

Use:
- exact existing P93 fastdiag message fingerprint;
- same encoded message/ciphertext construction;
- current Fast evaluator;
- Standard evaluator with >=3 independent key/evaluation-key trials;
- decoded numerical metrics over all slots.

No performance optimization.
No production arithmetic changes.
Do not reduce this to a simple 1e-2 pass/fail; report RMSE, quantiles, max/worst slot, precision bits, and Standard trial variability.

Return the required numerical-reference classification/metrics and `READY_FOR_WEB_REVIEW`.
