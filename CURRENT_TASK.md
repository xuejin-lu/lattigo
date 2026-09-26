# Current Task

Task: QPREFIX-AUDIT-002
Status: READY_FOR_CODEX

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap@main`
`specs/QPREFIX-AUDIT-002-C2S-RAW-RESCALE.md`

Task class:
`M/E — Focused evidence completion`

Audit current `fast-qprefix` C2S only.

If necessary, add a focused diagnostic `*_test.go` that manually steps:
`LinearTransform -> Rescale -> restore`
for all four current LogN13/P93 C2S groups and records exact bounds/levels/scales.

Production source must remain unchanged.
Do not modify `fast-ckks`, EvalMod, S2C, parameters, or schedules.

Commit/push any useful diagnostic test only after validation.
