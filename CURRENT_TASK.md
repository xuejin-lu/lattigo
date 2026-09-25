# Current Task

Task: FAST-INTEGRATION-001
Status: READY_FOR_CODEX

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap@main`
`specs/FAST-INTEGRATION-001-PRIVATE-F-MODUP-BASIS-BRIDGE.md`

Primary task pointer:
`heart-lattigo-bootstrap/CURRENT_TASK.md`

Task class:
`I — Implementation`

Accepted prerequisites:
- FAST-STORAGE-005 at `531aca50b5b38741e4e71cc98ea4b626bf88cb84`.
- Fixed-width-3 production policy at `d9919f9c080e0dfa731746f5c447f93633ae2f36`.

Implement only the first bounded production integration seam defined by the Primary spec:
- replace the historical basis-raise portion of `FastEvaluator.modUpBasis` with Level-0 LogicalQ -> ImportLevel0(width=3) -> FastStorageModUpLevel0(MaxLevel) -> compact maintained LogicalQ;
- add the compact private-F -> logical-Q bridge that materializes only existing maintained q rows;
- preserve existing scale alignment, Trace, Montgomery conversion, downstream DFT/EvalMod/S2C, packing, and production Rescale semantics;
- record ModUp-basis benchmark evidence.

Do not migrate downstream arithmetic to `FastCiphertext`.
Do not add KeySwitch, Relinearize, Rotate, contraction/adaptive width, frontend flags, or full-RNS fallback.

Follow the normal bounded implementation -> self-review -> at most one repair -> validation workflow.
Commit and push Secondary `fast-ckks`, then report `READY_FOR_WEB_REVIEW`.
