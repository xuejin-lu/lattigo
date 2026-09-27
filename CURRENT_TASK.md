# Current Task

Task: QPREFIX-IMPL-004-REVIEW-FIX
Status: READY_FOR_CODEX

Accepted candidate:
`18f4da53f03acd065d18550a8f2106725462896f`

Revised Primary spec:
`specs/QPREFIX-IMPL-004-RESCALE-LEVEL-TRANSITIONS.md`
at
`a7bf5f2c68662f66c7b955b75a566e5eff4e0816`

Review blocker:
`DropLevelCanonical` must not accept targetLevel > 3.

Reason:
for targetLevel > 3, q0123 does not determine the canonical centered representative
modulo full logical Q_target.

Required bounded fix:
- reject canonical targetLevel > 3 before mutating output;
- add transactional 5->4 (or equivalent) rejection test;
- keep existing <=3 canonical behavior unchanged;
- do not modify SameLift or any other subsystem.

Run focused and full regressions, commit/push, then report `READY_FOR_WEB_REVIEW`.
