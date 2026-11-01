# Fast CKKS public-key zero-secret API compatibility

User-authorized implementation scope: Public Key API Compatibility (方案 A).
Target: `experiment/fast-bootstrap-noise-fidelity`.
Parent: `be9c42ff23662300769b77f1eead3f12fc53adef`.

## Contract and implementation

`core/rlwe/encryptor.go` dispatches non-nil plaintext `Encrypt` / `EncryptNew`
with either a validated `*SecretKey` or `*PublicKey` to the existing
`encryptFastCKKSZeroSecret`, only when the parameter provider exposes
`FastCKKSZeroSecretSimulation` and the output is `*Ciphertext`.
The helper's degree-one storage, ring degree, metadata and level checks are
unchanged. It copies encoded plaintext into c0, clears c1, and preserves
Level/Scale/NTT/Montgomery metadata with independent output ownership.
This is intentionally insecure numerical simulation; native PK randomness and
noise are not preserved. Key generation, EncryptZero, nil-plaintext Encrypt,
and evaluation-key/invariant checks are unchanged.

## Validation

Go 1.26.4, darwin/arm64. Commands use
`GOCACHE=/private/tmp/lattigo-pk-bridge-go-cache` because the default cache is
not writable in the execution sandbox.

Passed:

- `go test ./schemes/ckks -run '^(TestFastCKKSZeroSecret|TestCKKSPublicKeyEncryptionWithoutFastCapability|TestGenericFastCKKSKeyGenerator|TestFastQPCompatibility)' -count=1`
- `go test ./schemes/ckks -run '^TestCKKS$/Evaluator/MulRelinThenAdd' -short -count=1`
- `go test ./core/rlwe ./schemes/bgv -short -count=1`
- `go vet ./core/rlwe ./schemes/ckks ./schemes/bgv`
- `git diff --check`

Coverage includes SK and PK EncryptNew, exact c0/plaintext equality and zero
c1 at every active row, decoded complex values (1e-7 tolerance), Levels 0/1,
coefficient/NTT domains, sparse/full slots, deep metadata/Scale ownership,
preallocated Encrypt at a lower output level, malformed output rejection,
WithKey, native EncryptZero and nil-plaintext Encrypt, no-key rejection, and
existing invalid-PK-layout rejection. PK ciphertexts pass ordinary RotateNew
and the exact DecomposeNTT / AutomorphismHoistedLazy / ModDownQPtoQNTT path
with numerical oracles and Fast key-switch-elision counters.

The new CKKS provider-isolation regression hides the Fast capability and
checks that native PK c1 remains nonzero and decoding succeeds. This test in
the Fast checkout is not labeled a genuine Standard baseline.
Separately, original Standard source was exported from
`5dbffbdea05394de2ca3a432ed5318aa832e3f40` into an isolated temporary directory;
only the new isolation test was added, with `nativeParams := params` so that
ordinary Standard CKKS parameters are passed directly to the public APIs.
`GOWORK=off go test ./schemes/ckks -run '^TestCKKSPublicKeyEncryptionWithoutFastCapability$' -count=1`
passed there without modifying Standard production source.

Before the fix, the new PK plaintext-preservation test failed on the original
SK-only dispatch; after the fix, both SK and PK cases pass. The existing
MulRelinThenAdd rejection test now explicitly injects invalid nonzero c1,
instead of assuming that PK encryption necessarily creates invalid Fast input.
Its transactional rejection and evaluator guard remain covered.

## Existing limitations observed during broader checks

The whole CKKS suite is not green. These unrelated areas were not repaired:

- `TestPublicFastZeroSecretAddFallsBackForFullCiphertextAndRejectsCompactGenericOperand`
  expects an error containing `missing active row q4` for compact scalar Add,
  but gets nil. The exact test also fails in an isolated export of parent
  `be9c42ff...`, at the same assertion (line 585).
- `go test ./schemes/ckks -short -count=1` fails in both the parent export and
  this checkout. Both show stale degree-0 error-text expectations and
  RescaleTo/Many precision failures. Parallel tests also encounter unsupported
  decrypt-projection panics: the parent can stop on nonzero non-c0 components;
  this checkout can stop on the unsupported ConjugateInvariant compact ring.
  This task does not establish whole-suite or ConjugateInvariant acceptance.
- `RotateHoistedNew` rejects Fast-layout Galois keys through
  `AutomorphismHoisted` (`cannot use Fast evaluation key with Standard execution`).
  The failure was also reproduced with SK-encrypted input, whose encryption
  path is unchanged. This separate API limitation is not the reported CNN512
  `AutomorphismHoistedLazy` / nonzero-c1 failure. The latter path passes.

No CNN512 integration, Bootstrap call, performance benchmark, Primary change,
key-generator change, or noise-model change was performed. The library is
ready for behemoth1 to update this branch and recompile CNN512 for its remote
integration test; this is not a claim that the full CNN has already passed.

Coding self-review checked the complete diff, capability/key/output guards,
helper reuse, metadata ownership, native-path isolation and retained c1
rejection. The only compatibility-test repair was the stale implicit PK-c1
assumption above. Scientific acceptance remains `READY_FOR_WEB_REVIEW`.
