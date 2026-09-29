# Current Task

Task: QPREFIX-PERF-OPT-003
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/QPREFIX-PERF-OPT-003-SOURCE-INTT-BATCHING.md`

Task class:
`P — Performance Repair`

Purpose:
Improve source-domain INTT locality by transforming active rows in row-major order across ciphertext components, using existing evaluator-owned staging.

Hard boundaries:
- no Q-prefix width/policy changes;
- no arithmetic/rounding/capacity changes;
- no transactional weakening;
- no goroutines, unsafe, assembly, or ring-kernel rewrite;
- no P93/generated-power schedule changes;
- no F/full-RNS fallback.

Run feasibility first. If rows4 source conversion improves <3%, do not retain production churn; report no-win.

Return the required optimization classification and `READY_FOR_WEB_REVIEW`.
