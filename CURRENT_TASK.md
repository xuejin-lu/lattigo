# Current Task

Task: QPREFIX-IMPL-003
Status: READY_FOR_CODEX

Authoritative specification:
Primary `specs/QPREFIX-IMPL-003-PREFIX-COMPLETE-PRIMITIVE-KERNELS.md`

Accepted prerequisites:
- `6553491f9fb9b964c8fd0d743f3a302de83d0b54`
- `91baa6a4655e10fe2460a399633fa54a03a35318`

Implement prefix-complete arithmetic/domain kernels for explicit widths 1..4.

Do not equate allocated backing with authoritative rows.
Do not globally replace legacy maintained-row policy in production wrappers yet.
Poisoned-row transition tests are mandatory.

Do not touch Rescale, ModUp, DFT/EvalMod/Bootstrap semantics or `fast-ckks`.

Commit/push `fast-qprefix` and report `READY_FOR_WEB_REVIEW`.
