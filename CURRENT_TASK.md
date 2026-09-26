# Current Task

Task: QPREFIX-AUDIT-001
Status: READY_FOR_CODEX

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap@main`
`specs/QPREFIX-AUDIT-001-FULL-BOOTSTRAP-CAPACITY.md`

Task class:
`M/E — Architecture audit + reproducible measurement`

This branch intentionally starts from pre-F commit:
`40532b4dce5c7eeae2db5b0b6f21be64801ce923`

Authoritative architecture:
`docs/FAST_QPREFIX_SPEC.md`

Audit only:
- no production behavior changes;
- no independent F basis;
- maintained prefix policy is `q0..q[min(level,3)]`;
- prove or disprove strict centered capacity at every required Bootstrap stage;
- identify actual Level-crossing contraction/canonicalization semantics.

Use test/diagnostic instrumentation only if required and keep it isolated.

Do not modify `fast-ckks`.
