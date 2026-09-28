# Current Task

Task: QPREFIX-PERF-DIAG-001
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-PERF-DIAG-001-ATTRIBUTION.md`

Production candidate:
`f9c7f21e65915bd3eafcd5b12590b570c22a7d6f`

Historical baseline:
`40532b4dce5c7eeae2db5b0b6f21be64801ce923`

This is diagnosis only.

Important:
- historical matched q0=55 baseline generally runs legacy q01 authority at high Level;
- current Q-prefix v2 production runs exact `QPrefixWidth(Level)`, i.e. q0123 at Level>=3;
- do not call the ~3.8x result a simple q012->q0123 implementation regression.

Measure and attribute:
- row-width scaling 2/3/4;
- Rescale CRT/division;
- polynomial/PS phases;
- one-bit guard;
- DoubleAngle;
- copy/allocation traffic;
- DFT/C2S/S2C;
- packing/ring-degree.

Do not modify production code.
Do not alter schedules/parameters/authority policy.
Do not introduce F/full-RNS fallback.

Temporary benchmark/test files are allowed only for measurement and must be removed before completion.

Return the required diagnostic status and `READY_FOR_WEB_REVIEW`.
