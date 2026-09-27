# Current Task

Task: QPREFIX-IMPL-008
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-IMPL-008-PUBLIC-BOOTSTRAP-BOUNDARY.md`

Accepted prerequisite:
`74c058ad59655f2a47efcb4faf1cf38324bd6137`

Implement only the public/structural boundary:
- explicit-row N1<->N2 conversion up to 4 rows;
- explicit-row packing/unpacking;
- production BootstrapMany public authority selection;
- finalization/public output contract;
- dormant-row isolation and public API equivalence.

Important:
- current public Residual MaxLevel <=1, so production public rows are only 1 or 2;
- q0123 remains an internal Bootstrap representation;
- do not widen public ciphertexts just because internal circuit is q0123;
- keep legacy ring-degree/packing wrappers legacy-authority.

Do not alter C2S/EvalMod/S2C mathematics, polynomial/DFT schedules, parameters, F/full-RNS fallback, or `fast-ckks`.

Commit/push `fast-qprefix` after focused and full regressions pass, then report `READY_FOR_WEB_REVIEW`.
