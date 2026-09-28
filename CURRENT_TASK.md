# Current Task

Task: QPREFIX-PERF-OPT-001
Status: READY_FOR_CODEX

Working branch: `fast-qprefix`
Required AGENTS revision: `FAST-QPREFIX-ACTIVE-001`

Startup assertion:
- after sync, verify HEAD is on `fast-qprefix`;
- re-read `AGENTS.md`;
- confirm it contains `FAST-QPREFIX-ACTIVE-001`;
- if the running agent still believes `AGENTS.md` requires `fast-ckks`, treat the agent session as stale and stop with `STALE_AGENT_INSTRUCTIONS` rather than switching branches.

Authoritative Primary specification:
`specs/QPREFIX-PERF-OPT-001-TRANSACTIONAL-RESCALE-STAGING.md`

Task class:
`P — Performance Repair`

Purpose:
Optimize shared Q-prefix Rescale by replacing the current two-exact-computation transactional design with evaluator-owned staging:

1. compute exact target residues once for all components;
2. validate all capacity constraints before caller-visible mutation;
3. only then commit staged residues through NTT/Montgomery restore.

Do not:
- weaken failure-before-mutation;
- change CRT/rounding semantics;
- change Q-prefix width/policy;
- change P93 or polynomial schedule;
- introduce F/full-RNS fallback.

Use the reusable fastdiag hooks for post-change attribution and normal diagnostics-OFF benchmarks for authoritative speed.

Return the required optimization classification and `READY_FOR_WEB_REVIEW`.
