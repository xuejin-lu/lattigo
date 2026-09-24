# Current Task

Task: FAST-STORAGE-002
Status: AUTHORIZED_BY_PRIMARY

Primary specification:
`xuejin-lu/heart-lattigo-bootstrap/specs/FAST-STORAGE-002-CONTAINER-CONVERSION-BOUNDARIES.md`

Scope:
Implement and test the physically separate Fast storage ciphertext container and Level-0 LogicalQ <-> FastStorage conversion boundaries.

Architecture authority:
`docs/FAST_CKKS_SPEC.md`

Important:
- f_i rows must never masquerade as logical q_i rows;
- LogicalLevel is explicit and independent from ActiveStorageWidth;
- do not wire this container into Bootstrap, Rescale, ModUp, arithmetic, or keys in this task.

Follow `AGENTS.md`, safely sync `fast-ckks`, read the Primary spec, then execute only the authorized scope.
