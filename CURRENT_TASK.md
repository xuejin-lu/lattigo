# Current Task

Task: QPREFIX-IMPL-005
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-IMPL-005-LEVEL0-MODUP-TRACE.md`

Accepted prerequisite:
`18f4da53f03acd065d18550a8f2106725462896f`

Implement only the Level-0 Bootstrap entry migration:
- canonical q0 ModUp with positive midpoint;
- exact q-prefix materialization through q3 when Level >=3;
- explicit-width integer scale alignment;
- explicit-width Trace;
- consistent final Montgomery conversion across all authoritative rows.

ScaleDown stays legacy-authority before ModUp.
DFT/C2S remains legacy q012 and may intentionally drop q3 authority after this boundary until QPREFIX-IMPL-006.

Do not expand into DFT/LinearTransform, EvalMod/PS/DA, packing/N1-N2, parameters, F, full-RNS fallback, or `fast-ckks`.

Commit/push `fast-qprefix` and report `READY_FOR_WEB_REVIEW`.
