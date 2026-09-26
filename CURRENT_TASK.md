# Current Task

Task: QPREFIX-AUDIT-003
Status: READY_FOR_CODEX

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap@main`
`specs/QPREFIX-AUDIT-003-EVALMOD-PS-DA-CAPACITY.md`

Task class:
`M/E — Focused evidence completion`

Audit current `fast-qprefix` EvalMod/PS/DoubleAngle only.

If necessary, add diagnostic/test-only instrumentation to expose:
- generated-power operands/products/rescales;
- PS baby/giant boundaries;
- DoubleAngle rounds;
- exact maintained-prefix coefficient bounds.

Production semantics must remain unchanged.
Do not modify `fast-ckks`, S2C, parameters, planScale, or polynomial scheduling.

Commit/push useful diagnostic-only instrumentation only after validation.
