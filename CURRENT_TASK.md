# Current Task

Task: QPREFIX-IMPL-009
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/QPREFIX-IMPL-009-PERFORMANCE-RELEASE-GATE.md`
at Primary commit:
`3ace6637ac793db1e246a1abfd9ccd2612abab1b`

Production candidate:
`f9c7f21e65915bd3eafcd5b12590b570c22a7d6f`

008 constitutional re-review:
production authority is NOT "internal always q0123".

The exact production invariant is:
`rows = QPrefixWidth(Level) = min(Level+1,4)`.

Therefore:
- L0 -> q0
- L1 -> q01
- L2 -> q012
- L>=3 -> q0123
- after every logical Level change, authority must be recomputed/contracted accordingly.

Low-level explicit-row helpers may support narrower widths for legacy/transition compatibility, but production Q-prefix v2 must use the exact Level-defined width.

The release gate must record the full production stage map and assert equality, not merely `rows <= 4`.

Continue QPREFIX-IMPL-009 validation only.
Do not change production source unless a hard release bug is found and reported first.
