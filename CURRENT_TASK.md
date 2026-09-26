# Current Task

Task: QPREFIX-IMPL-002
Status: READY_FOR_CODEX

Authoritative architecture:
- `docs/FAST_QPREFIX_SPEC.md`
- Primary `docs/QPREFIX-V2-PRODUCTION-MIGRATION-PLAN.md`

Accepted prerequisite:
- QPREFIX-IMPL-001 at `6553491f9fb9b964c8fd0d743f3a302de83d0b54`.

Task class:
`I — Implementation`

Implement only compact ciphertext lifecycle integration:
- constructors;
- copy/copy-new;
- physical resize/backing width;
- degree/output allocation;
- evaluator-owned scratch.

Use the single `QPrefixWidth(level)` policy.

Do not widen arithmetic producers/consumers yet.
Do not change representative semantics in structural Resize.
Do not touch Rescale, ModUp, DFT, EvalMod, Bootstrap, `fast-ckks`, or add F/adaptive width.

Required validation:
- focused lifecycle tests;
- relevant fast package tests;
- `go test ./...`;
- `git diff --check`.

Commit/push `fast-qprefix` and report `READY_FOR_WEB_REVIEW`.
