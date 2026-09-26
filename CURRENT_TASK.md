# Current Task

Task: FAST-STORAGE-008
Status: READY_FOR_CODEX

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap@main`
`specs/FAST-STORAGE-008-C2S-WIDTH2-EXPERIMENT.md`

Task class:
`E — Experiment / bounded implementation`

Accepted prerequisites:
- FAST-STORAGE-007 at `1a8018efda159f1dc9ab7078ce80a40c6db28e71`.
- Explicit stage-local width experiment rule at `04396c94a3bb063e93b00bb80a500dfcc2891e21`.

Implement only:
- explicit storage contraction 3->2 under strict bound proof;
- width-2 plaintext mirror and LinearTransform support;
- actual LogN13 C2S capacity experiment;
- exact Logical/F3/F2 comparisons;
- representative factor and full-chain benchmarks.

Do not modify production Bootstrap/C2S or enable adaptive runtime width.

Commit/push `fast-ckks` and report `READY_FOR_WEB_REVIEW`, or `NEEDS_WEB_REVIEW` only for architecture/math conflicts.
