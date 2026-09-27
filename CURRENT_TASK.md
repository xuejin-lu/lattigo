# Current Task

Task: QPREFIX-IMPL-006
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-IMPL-006-LINEARTRANSFORM-DFT-C2S-S2C.md`

Accepted prerequisite:
`08b36b0b730dcc57eb98466594796c10bffbfdb9`

Implement:
- explicit-row LinearTransform direct + BSGS;
- explicit-row DFT factor execution;
- production C2S q0123 activation from the migrated ModUp boundary;
- exact q3 group/checkpoint/capacity evidence;
- S2C q-prefix capability without premature production q0123 activation.

Do not migrate EvalMod/PS/DoubleAngle.
Do not assume allocated q3 is authoritative at EvalMod->S2C.
Do not touch packing/N1-N2, parameters, F, full-RNS fallback, or `fast-ckks`.

Commit/push `fast-qprefix` after focused and full regressions pass, then report `READY_FOR_WEB_REVIEW`.
