# Current Task

Task: DIAG-FRAMEWORK-001
Status: READY_FOR_CODEX

Authoritative Primary specification:
`specs/DIAG-FRAMEWORK-001-REUSABLE-TRACE.md`

Task class:
`I — Diagnostic Infrastructure`

Purpose:
Add compile-time-gated reusable Fast diagnostic hooks for exactly these trace scopes:

- `stage`
- `power`
- `rescale`

Normal builds must compile diagnostics away and must not add persistent hot-path timing/allocation overhead.

The framework must preserve arithmetic and transactional semantics.

Do not:
- optimize or fuse Rescale;
- remove preflight;
- change Q-prefix arithmetic;
- change `QPrefixWidth(Level)`;
- change generated-power/P93 scheduling;
- introduce F/full-RNS fallback.

Follow the Primary spec for exact event schema, build-tag behavior, tests, and completion criteria.

Return the required framework classification and `READY_FOR_WEB_REVIEW`.
