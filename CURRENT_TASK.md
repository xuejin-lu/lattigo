# Current Task

Task: FAST-INTEGRATION-003
Status: READY_FOR_CODEX

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap@main`
`specs/FAST-INTEGRATION-003-PRIVATE-F-RESIDENT-MODUP-TRACE.md`

Task class:
`I — Implementation`

Accepted prerequisites:
- FAST-INTEGRATION-002 at `57ffb88744c82778c0a9392ecab394e19f712a3d`.
- Normalized private-F Trace theorem at `8f823fdb9464c2738f71c30d156ce574098d8605`.

Implement only the resident private-F ModUp + scale-alignment + normalized-Trace segment defined by Primary.

Key frozen rules:
- production private width = 3;
- Trace computes unnormalized automorphism sum first, then exact normalization;
- require `2*g*B < S3` for every component before Trace;
- no private-F pre-multiplication by modular `g^-1`;
- no LogicalQ/full-RNS fallback;
- export compact LogicalQ only after Trace.

Preserve downstream Montgomery/DFT/EvalMod/packing/production Rescale.

Commit and push Secondary `fast-ckks`, then report `READY_FOR_WEB_REVIEW` or `NEEDS_WEB_REVIEW`.
