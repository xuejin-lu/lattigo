# Current Task

Task: QPREFIX-IMPL-007-REVIEW-FIX
Status: READY_FOR_CODEX

Accepted candidate:
`980be10e7ad8d14df4541526c05b4bcb44f17996`

Only blocker:
legacy exported polynomial/Mod1 wrappers must not infer q0123 authority from allocated backing.

Required bounded fix:
- restore legacy row authority on existing exported polynomial Evaluate wrappers;
- add explicit-row polynomial entry points and keep all current q0123 implementation behind them;
- restore legacy authority on `mod1.FastEvaluator.EvaluateNew`;
- add explicit-row Mod1 entry point;
- production Bootstrap EvalMod explicitly passes `QPrefixWidth(input.Level())` to that new path;
- keep current q0123 capacity/PS/guard/DoubleAngle/S2C behavior unchanged.

Add q3 poison compatibility tests proving:
- legacy wrappers ignore non-authoritative q3;
- explicit rows=4 paths consume q3.

Do not modify the polynomial schedule, parameters, capacity logic, S2C mathematics, packing/N1-N2, or `fast-ckks`.

Run focused/full regressions and push normally, then report `READY_FOR_WEB_REVIEW`.
