# Current Task

Task: FAST-INTEGRATION-002
Status: READY_FOR_CODEX

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap@main`
`specs/FAST-INTEGRATION-002-FUSED-MODUP-BRIDGE-OPTIMIZATION.md`

Task class:
`I — Implementation`

Accepted prerequisite:
- FAST-INTEGRATION-001 at `31efadc693559217b48d3e76a2e3655b9e6cd14d`.
- Fusion architecture clarification at `dc698e2d99a488f0de2cf4f3207b09ef94521303`.

Implement only the fused Level-0 ModUp bridge optimization defined by Primary.

Preserve accepted canonical q0 semantics, compact maintained logical rows, fixed-width-3 architecture, standalone private-F reference APIs, and downstream Bootstrap behavior.

Performance gates are part of task acceptance.

Commit and push Secondary `fast-ckks`, then report `READY_FOR_WEB_REVIEW` or `NEEDS_WEB_REVIEW`.
