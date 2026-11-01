# Fast Experimental Trust Mode — CNN512 compatibility

User decision: 2026-11-02. Branch: `experiment/fast-bootstrap-noise-fidelity`.
Parent: `efba6cb594d6c09b2d93e08db32439cdfd248383`.
Status: `CAPACITY_UNPROVEN_EXPERIMENTAL`, `READY_FOR_WEB_REVIEW`.

## Policy and reuse

This branch defaults ordinary Fast-capability CKKS evaluators to Experimental
Trust for the existing QP compatibility lane. There is no new evaluator,
application flag, environment selector, or CNN source modification. Original
Standard and generic RLWE do not enter this capability-gated lane.

A fresh evaluator can perform high-level ModDown with complete caller-supplied
Q/P backing without inheriting `armed`. The existing basis extender performs
the actual ModDown. The Full-Q helpers still validate backing; Trust skips
fingerprint hashing, INTT/CRT centered-prefix consistency, P-inverse validation
and extra capacity proofs. It does not mark `armed` or certify a ciphertext.
The immutable internal policy can be disabled in strict regression fixtures;
the old strict checks and numerical thresholds are retained.

Ciphertext/QP metadata, dimensions, Q/P levels and required row backing remain
checked at their existing API boundaries. Zero-secret prefix checks, zero P
for ModDown, Fast-key shape/layout, decomposition checks and pending-error
propagation remain active. Missing input rows are rejected before full-RNS
operations; no input rows are invented. Zero-key-switch shortcuts are unchanged.

Complete high-level, all-zero-c1 ciphertexts use the existing full-Q Rescale
kernel even on a new evaluator. Strict policy retains the historical armed
gate. Compact ciphertexts and prefix-owned inputs with dormant nonzero c1
rows keep the existing shared Q-prefix Rescale core. Complete storage alone
is not a numerical certificate: supplied c0 rows can be incoherent, and Trust
accepts this explicitly authorized numerical risk. No blanket full-RNS fallback
was added to other Fast arithmetic or to the explicit Fast Bootstrap engine.

## Lightweight validation

Go 1.26.4, darwin/arm64; `GOCACHE=/private/tmp/lattigo-pk-bridge-go-cache`.

The new test uses LogN10, MaxLevel26, five P primes, L15 / scale 2^32, and the
CNN chain's 42/32-bit shape. Smaller N keeps this regression lightweight; it is
not the full LogN16 CNN integration. Evaluator A performs baby rotation 1;
two actual Montgomery mask multiply/add terms accumulate in full QP. A fresh
Evaluator B performs c1 ModDown, the zero gadget shortcut, giant rotation 64,
and final c0/c1 ModDown. A third fresh evaluator performs actual full-Q Rescale.
The same fixture also exercises shared compact Q-prefix Rescale.

Observed maximum complex errors against the plaintext oracle (limit 1e-6):

| Checkpoint | Maximum error |
|---|---:|
| L15 accumulated ModDown, 0.5 * Rotate(input, 65) | 1.71617286025e-9 |
| L14 full-Q Rescale | 5.20151746355e-9 |
| L14 shared compact Q-prefix Rescale | 5.20151746355e-9 |

Full-Q Rescale matches the independent full-RNS division oracle on every row;
Scale/Level, zero c1, unchanged source, distinct unarmed states, empty trusted
fingerprint cache, key-switch/gadget elisions and compact routing are asserted.
The strict-policy fresh evaluator still rejects the same unarmed L15 ModDown.
Negative cases retain missing-Q/P/output-row, invalid-Level, nonzero-P and
nonzero authoritative gadget-input rejection without output mutation. A test
separately proves that incoherent high rows are accepted by Trust, rejected by
Strict, while missing high rows remain rejected by both.

Passed:

- 34 top-level CKKS tests selected by
  `go test ./schemes/ckks -run '^(TestFastQP|TestFastCKKSZeroSecret|TestGenericFastCKKSKeyGenerator|TestCKKSPublicKeyEncryptionWithoutFastCapability|TestPublicFastAutomorphism|TestPublicFastBTS|TestPublicFastZeroSecretCompactAddMulRescaleRotateChain|TestPublicFastZeroSecretRotate)' -count=1 -v`.
- `go test ./core/rlwe ./schemes/bgv -short -count=1`.
- `go test ./circuits/ckks/bootstrapping -run '^(TestFastE32GenEvaluationKeysAndPublicEvaluatorDispatch|TestNewEvaluatorFastDispatchDoesNotRequireDenseSparseKeys)$' -count=1` (constructor/dispatch only; zero Bootstrap calls).
- `go vet ./schemes/ckks ./core/rlwe ./schemes/bgv` and `git diff --check`.
- Genuine original Standard `5dbffbdea05394de2ca3a432ed5318aa832e3f40`, exported
  without production changes to an isolated temporary directory:
  `GOWORK=off go test ./schemes/ckks ./core/rlwe -short -count=1`.

A broader CKKS selection still reproduces the known unrelated
`TestPublicFastZeroSecretAddFallsBackForFullCiphertextAndRejectsCompactGenericOperand`
failure at line 585 (expected compact scalar Add error, got nil). Its parent
reproduction and other whole-suite debt are already recorded in
`FAST_PUBLIC_KEY_ZERO_SECRET_COMPATIBILITY.md`; no whole-suite acceptance is
claimed and that test was not changed or suppressed.

Coding self-review checked the complete diff, strict-policy preservation,
capability isolation, actual arithmetic, full/compact Rescale selection and
structural-error propagation. The negative-test snapshot was repaired to
preserve nil row identity because Poly.CopyNew normalizes nil to empty slices.
No numerical failure or threshold relaxation occurred. No benchmark, CNN,
Bootstrap, Primary, model, weight, Golden or noise-hook modification occurred.

behemoth1 can fetch this branch and recompile unchanged CNN512. Its complete
LogN16 Golden comparison, accumulated coefficient ranges and end-to-end
performance remain unmeasured here. The pre-existing WithKey helper/state-copy
limitation is outside this task; CNN's observed path uses NewEvaluator.
