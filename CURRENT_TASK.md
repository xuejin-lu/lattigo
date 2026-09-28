# Current Task

Task: DIAG-FRAMEWORK-001-R1
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/DIAG-FRAMEWORK-001-R1-RESCALE-EVENT-CLOSURE.md`

Task class:
`I — Diagnostic Infrastructure Repair`

Purpose:
Repair Rescale trace event nesting so diagnostic closure does not double-count the same reconstruction interval.

Required shape:
- Rescale parent
- preflight/materialization children
- component-level prefix/coefficient-loop/restore children
- reconstruction and residue events nested under coefficient-loop

Do not:
- change Rescale arithmetic;
- remove/fuse preflight;
- alter transactional semantics;
- optimize production code;
- change Q-prefix width/policy;
- change generated-power schedule;
- introduce F/full-RNS fallback.

Return the required repair classification and `READY_FOR_WEB_REVIEW`.
