# Current Task

Task: QPREFIX-IMPL-004
Status: READY_FOR_CODEX

Authoritative specification:
Primary `specs/QPREFIX-IMPL-004-RESCALE-LEVEL-TRANSITIONS.md`

Accepted prerequisites:
- `6553491f9fb9b964c8fd0d743f3a302de83d0b54`
- `91baa6a4655e10fe2460a399633fa54a03a35318`
- `c8b591a30c05a2261de8d0181d7b3f64bc58169b`

Implement only:
- q0123 fixed-width centered reconstruction;
- Rescale with source width from QPrefixWidth(level);
- logical q_level divisor even when level > 3;
- target prefix shrink with new Level;
- explicit SameLift and Canonical DropLevel semantics;
- strict transactional capacity failures.

Do not globally activate q3 outside Rescale/DropLevel.
Do not modify ModUp/Trace, DFT/EvalMod/Bootstrap, `fast-ckks`, parameters, or add F.

Commit/push `fast-qprefix` and report `READY_FOR_WEB_REVIEW`.
