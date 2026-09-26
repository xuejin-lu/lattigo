# Current Task

Task: QPREFIX-ARCH-PLAN-001
Status: READY_FOR_CODEX

Task class:
`M — Architecture planning only`

This branch is the architecture-reading target only.

Authoritative architecture:
`docs/FAST_QPREFIX_SPEC.md`

Important decision:
- Q-prefix v2 is the production direction;
- private F is not the production direction;
- `QPREFIX-AUDIT-003` is paused and must not be resumed unless a later Primary task explicitly re-authorizes it.

The user will provide the detailed architecture-planner prompt directly in Codex chat.

Do not modify this Secondary branch.
Do not write implementation code.
Do not run a broad repository-wide implementation audit beyond what the user's planning prompt requires.

The only repository write for this task belongs in Primary:
`docs/QPREFIX-V2-PRODUCTION-MIGRATION-PLAN.md`
