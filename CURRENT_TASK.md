# Current Task

Task: QPREFIX-IMPL-007
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-IMPL-007-EVALMOD-PS-DOUBLEANGLE.md`

Accepted prerequisite:
`2aaec605951b4469cef6db10c28453422204595e`

Implement:
- explicit q0123 authority through Mod1 input clone;
- polynomial workspace/power generation;
- PS baby/giant steps and scale alignment;
- explicit-row Mul/MulRelin/Relinearize/MulThenAdd/scalar Add/Sub;
- q0123 one-bit guard;
- all EvalMod Rescale boundaries;
- DoubleAngle;
- final EvalMod output q0123 contract;
- production S2C activation only after that contract is proven.

For accepted P93, Level remains >=4 through EvalMod, so rows=4 must remain authoritative throughout.

Preserve the existing accepted P93 polynomial/PS/scale/DoubleAngle schedule.
Do not redesign polynomial coefficients, PS split, parameters, DFT factorization, packing/N1-N2, or public error gates.
Do not add F/full-RNS fallback or modify `fast-ckks`.

Commit/push `fast-qprefix` after focused and full regressions pass, then report `READY_FOR_WEB_REVIEW`.
