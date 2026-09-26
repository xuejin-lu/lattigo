# Current Task

Task: FAST-STORAGE-007
Status: READY_FOR_CODEX

Authoritative Primary specification:
`xuejin-lu/heart-lattigo-bootstrap@main`
`specs/FAST-STORAGE-007-PRIVATE-F-PLAINTEXT-LINEAR-TRANSFORM.md`

Task class:
`I — Implementation / feasibility foundation`

Accepted prerequisites:
- FAST-STORAGE-006 at `532319346d8235fb42c72bfd22b57a6468675c82`.
- Private-F plaintext mirror architecture at `22b9f07969af38705573686fe96a873e4bd001a3`.

Implement only:
- private-F plaintext mirror from complete logical-Q encoded diagonals;
- exact plaintext bounds/L1 norms;
- private-F plaintext multiply and automorphism;
- standalone private-F LinearTransform;
- actual LogN13 C2S capacity audit.

Do not modify production Bootstrap/C2S or matrix generation.

If actual C2S width-3 capacity fails, stop with `NEEDS_WEB_REVIEW` and report the first failing factor + exact bound numbers.
If feasible, commit/push Secondary `fast-ckks` and report `READY_FOR_WEB_REVIEW`.
