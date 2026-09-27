# Current Task

Task: QPREFIX-IMPL-004
Status: READY_FOR_CODEX

Authoritative revised Primary specification:
`specs/QPREFIX-IMPL-004-RESCALE-LEVEL-TRANSITIONS.md`
at Primary commit:
`87125e13570a4e88dc26343d88810c23966768f7`

Important correction:
q0123 Rescale **capability** is implemented now, but production activation remains legacy-width until upstream producers are migrated.

Current C2S LinearTransform still owns q012 only. Therefore:
- production Rescale/RescaleTo must continue selecting legacy authoritative rows;
- explicit-width kernel/tests must support q0123;
- q3 backing must never be treated as authority automatically.

Controlled recovery:
local uncommitted changes from the stopped first 004 attempt are expected.
Do not reset/clean/stash them away indiscriminately.
Continue from them and revise toward the updated contract.

Do not modify LinearTransform/DFT production routing, ModUp/Trace, EvalMod/PS/DA, Bootstrap orchestration, or `fast-ckks`.

Commit/push `fast-qprefix` after all focused and full regressions pass, then report `READY_FOR_WEB_REVIEW`.
