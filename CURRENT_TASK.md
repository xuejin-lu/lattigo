# Current Task

Task: QPREFIX-IMPL-004
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-IMPL-004-RESCALE-LEVEL-TRANSITIONS.md`

Accepted prerequisite:
QPREFIX-IMPL-003 at
`c8b591a30c05a2261de8d0181d7b3f64bc58169b`

Implement:
- fixed-width q0123 centered reconstruction;
- Q-prefix Rescale/RescaleTo;
- logical q_ell divisor even when ell > 3;
- 3->2 / 2->1 / 1->0 target-prefix contraction;
- explicit SameLift vs Canonical DropLevel semantics;
- transactional failures;
- required Standard/math-big oracles and C2S regression.

Only Rescale and explicit Level-transition APIs become Q-prefix-policy authoritative in this milestone.

Do not globally activate q0123 in other production subsystems.
Do not touch ModUp/Trace, DFT/LinearTransform routing, EvalMod/PS/DA, Bootstrap orchestration, `fast-ckks`, or add F/full-RNS fallback.

Commit/push `fast-qprefix` and report `READY_FOR_WEB_REVIEW`.
