# Current Task

Task: QPREFIX-IMPL-007-REVIEW-FIX-2
Status: READY_FOR_CODEX

Accepted candidate:
`4c6d7f207ddbaa531eb6d1010ffa6cd8ab08eace`

Final blocker:
the real `bootstrapCore` bypasses the q0123-aware Bootstrap SlotsToCoeffs wrapper and still calls legacy `DFTEvaluator.SlotsToCoeffsNew`.

Required fix only:
- at the actual EvalMod->S2C production handoff, pass explicit `QPrefixWidth(EvalModOutput.Level())`;
- use the q-prefix-capable S2C API from QPREFIX-IMPL-006;
- keep the legacy DFT `SlotsToCoeffsNew` wrapper unchanged;
- add a regression through the real Bootstrap/BootstrapMany call graph proving q3 is consumed by production S2C.

Do not change polynomial/PS/DoubleAngle mathematics, capacity logic, DFT factorization, packing/N1-N2, parameters, or `fast-ckks`.

Run focused/full regressions, push normally, then report `READY_FOR_WEB_REVIEW`.
