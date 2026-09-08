# Fast-CKKS repository guidance

## Purpose and source of truth

This `fast-ckks` branch is an intentionally insecure Fast-CKKS backend built on Lattigo. Keep application CKKS usage unchanged and keep Normal Lattigo behavior available. The current engineering priority is speed: Fast execution must become substantially faster than Normal Lattigo for large workloads.

Before Fast-CKKS work, read [`docs/FAST_CKKS_SPEC.md`](docs/FAST_CKKS_SPEC.md) and inspect the current source for the requested path. The specification is the single detailed source of truth; repository evidence takes precedence over API names or pretrained assumptions.

## Persistent invariants

- `q0` and `q1` are the authoritative RNS limbs in Fast hot paths. `q2...qL` are dormant/non-authoritative; do not silently perform full-RNS work or Standard/full-RNS fallback.
- `c0`, `c1`, `c2`, ... name ciphertext components. Do not confuse them with `q0`, `q1`, `q2`, ... RNS limbs.
- The current zero-secret implementation is a mode, not a permanent scientific invariant. Preserve extension points for sampled secrets and future noise experiments.
- Preserve public APIs and CKKS parameter objects where practical; preserve structure while eliding expensive dormant or security-only computation.
- Do not add expensive noise-fidelity work during Stage A. Normal and Fast implementations must coexist for correctness and performance comparison.
- Keep changes CKKS-specific where possible. Do not modify BFV, BGV, TFHE, or unrelated infrastructure merely for consistency.

For detailed architecture, current implementation status, hot-path rules, validation, and future experimental knobs, use `docs/FAST_CKKS_SPEC.md`.

## Workflow

Trace the exact source call path, identify the smallest modification surface, change one bounded behavior at a time, run targeted tests, and review the diff. Do not perform broad source audits or repository-wide refactors without explicit scope.
