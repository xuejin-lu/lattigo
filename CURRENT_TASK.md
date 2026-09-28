# Current Task

Task: QPREFIX-PERF-OPT-002
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Authoritative Primary specification:
`specs/QPREFIX-PERF-OPT-002-FIXED-WIDTH-DIVISION-DEDUP.md`

Task class:
`P — Performance Repair`

Purpose:
Optimize the current Q-prefix Rescale fixed-width hot loop with exact division deduplication and prepared invariant reconstruction constants.

Hard boundaries:
- no Q-prefix width/policy change;
- no rounding/CRT/capacity semantic change;
- no transactional weakening;
- no reciprocal approximation;
- no unsafe/assembly;
- no P93/generated-power schedule change;
- no F/full-RNS fallback.

Run the required correctness, focused benchmark, E2E, and low-overhead power-trace gates.

Return the required optimization classification and `READY_FOR_WEB_REVIEW`.
