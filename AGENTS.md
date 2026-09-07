# Fast-CKKS repository guidance

## Project purpose

This `fast-ckks` branch is an experimental Fast-CKKS implementation built on Lattigo. Modify the backend Lattigo implementation while keeping application code unchanged. Fast mode is intentionally not cryptographically secure.

## Source of truth

Before designing or implementing Fast-CKKS behavior:

1. Read [`docs/FAST_CKKS_SPEC.md`](docs/FAST_CKKS_SPEC.md).
2. Inspect the actual current Lattigo source involved in the requested operation.
3. Prefer repository-source evidence over API names or pretrained knowledge.

Do not infer implementation behavior from function names alone.

## Strict scope

Primary scope is CKKS, CKKS bootstrapping, CKKS evaluator behavior, and directly required RLWE primitives. BFV, BGV, TFHE, unrelated schemes, unrelated examples, and repository-wide refactors are out of scope unless a CKKS call path proves a shared change unavoidable. Never modify unrelated schemes merely for consistency.

## Minimal-change rule

Use the smallest possible modification surface. Before changing shared `core/rlwe` code, determine whether the change can remain CKKS-specific. If shared code must change, document why CKKS requires it and check other schemes for unintended impact. Do not redesign Lattigo globally or duplicate large implementations unless technically unavoidable.

## Fast-CKKS invariants

- Fast mode is intentionally insecure; do not silently fall back to secure execution.
- Preserve the existing CKKS parameter configuration. Do not reduce the configured modulus/level chain for speed.
- Relevant key material is intentionally zeroed according to the specification; zero keys are an intentional state, not an error.
- Fast execution is intended to operate with only the two degree-1 ciphertext components, named `c0` and `c1`.
- Reserve `q_i` for RNS modulus primes or modulus-chain levels. Do not confuse ciphertext components (`c0`, `c1`) with CKKS/RNS levels (`Q`, `q_i`, `LevelQ`). A 20-level configuration remains a 20-level configuration even when Fast execution retains two ciphertext components.

## Development workflow

For implementation tasks: inspect the exact operation, trace its source-level call path, identify the minimum modification surface, change one bounded behavior at a time, run targeted CKKS tests, and review the diff for unrelated changes. Do not begin with a repository-wide exploration or refactor.

## Testing

Prefer targeted tests for changed code. Do not run broad repository tests unless shared infrastructure changes justify them. Normal CKKS behavior must remain unchanged unintentionally, and Fast-specific behavior requires explicit tests.

## Git discipline

Work on the current Fast development branch unless instructed otherwise. Do not rewrite history. Keep commits small and single-purpose, and do not commit unrelated formatting or cleanup.
