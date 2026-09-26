# Current Task

Task: QPREFIX-IMPL-001
Status: READY_FOR_CODEX

Authoritative architecture:
- `docs/FAST_QPREFIX_SPEC.md`
- Primary `docs/QPREFIX-V2-PRODUCTION-MIGRATION-PLAN.md`

Task class:
`I — Implementation`

Implement only the Q-prefix v2 policy/capacity foundation.

Frozen policy:
[
w_Q(ell)=min(ell+1,4)
]

Required:
- single Level -> maintained-width policy;
- actual q-prefix product helper;
- per-component bound/capacity contract;
- strict `2B < S_Q` predicate;
- explicit transactional capacity error;
- focused tests including exact strict-boundary cases.

Do not connect the new policy to constructors, arithmetic, Rescale, ModUp, DFT, EvalMod, or Bootstrap yet.
Do not modify `fast-ckks`.
Do not introduce F or adaptive width.

Commit/push `fast-qprefix` and report `READY_FOR_WEB_REVIEW`.
