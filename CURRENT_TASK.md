# Current Task

Task: QPREFIX-PERF-OPT-004
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/QPREFIX-PERF-OPT-004-BARRETT-FIXED-WIDTH-REDUCTION.md`

Task class:
`P — Performance Repair`

Purpose:
Test exact Barrett-Horner fixed-width modular reduction in the Q-prefix Rescale hot path using existing Lattigo Barrett primitives.

Hard boundaries:
- feasibility first;
- exact arithmetic only;
- no approximate reciprocal / floating point;
- no unsafe / assembly / goroutines;
- no Q-prefix width/policy change;
- no CRT/rounding/capacity semantic change;
- no transactionality weakening;
- no P93/generated-power schedule change;
- no F/full-RNS fallback.

If feasibility fails, remove temporary candidate code and leave production unchanged.

Return the required optimization classification and `READY_FOR_WEB_REVIEW`.
