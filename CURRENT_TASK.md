# Current Task

Task: FAST-STANDARD-NUMERICAL-DIAG-002
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/FAST-STANDARD-NUMERICAL-DIAG-002-CURRENT-QPREFIX-LOCKSTEP.md`

Task class:
`D — Numerical Correctness Diagnosis`

Purpose:
On the current clean Q-prefix production source, identify where Fast first materially diverges numerically from genuine Standard for the canonical P93 q0=55 / 4096-slot workload.

Do:
- stage-by-stage semantic lockstep;
- genuine Standard decryption;
- current Fast semantic decode;
- exact scale/level/rows audit;
- internal EvalMod bisect only if EvalMod is first material;
- S2C amplification check.

Do not:
- modify production arithmetic;
- change plan scale;
- change q0/Q/P parameters;
- tune thresholds after seeing data;
- optimize performance.

Return the required numerical-bisect classification/checkpoints and `READY_FOR_WEB_REVIEW`.
