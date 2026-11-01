# Fast ciphertext + plaintext Add compatibility

User-authorized scope: CNN512 `stem.square` -> AddVecs, based on
`2608f4cc0063117077e176d730b1cf4bf23e74ec`, branch
`experiment/fast-bootstrap-noise-fidelity`.

## Implementation and boundaries

Ordinary `ckks.Evaluator.Add` adapts a `*rlwe.Plaintext` through a borrowed
(non-copying) degree-zero `rlwe.Ciphertext` view of its existing RLWE element.
For eligible equal-scale, matching-domain zero-secret operands, it dispatches
to the existing `addSubFastCKKSZeroSecret` and shared `fastcore.AddSubCore`.
The core already supports differing degrees: it adds plaintext to c0 and keeps
the ciphertext's higher components. No arithmetic kernel was added or changed.

Only min(Level+1, 4) authoritative rows are validated/read/computed. Output
Level is min(ciphertext, plaintext, output Level), Scale is preserved, and the
output is compact. Plaintext high rows are neither validated nor consumed,
whether fully backed, missing, or numerically poisoned. Input metadata, NTT,
Montgomery/batching state, prefix backing, canonical residues, equal scales
and zero non-c0 components retain their existing checks. Empty storage and nil
metadata are rejected before access. Ordinary in-place output is supported.

The pre-existing eligibility guard and generic fallback remain in place for
inputs outside this zero-secret contract; in particular, native complete-Q
nonzero-c1 ciphertexts retain their prior generic Add behavior. Compact inputs
with mismatched scales/domains are rejected before generic full-Q arithmetic.
No new scale-alignment algorithm, Sub/plaintext API, direct-vector Add,
Experimental Trust policy, key semantics or frontend mode switch was added.

## Validation evidence

Go 1.26.4, darwin/arm64. Commands used
`GOCACHE=/private/tmp/lattigo-pk-bridge-go-cache` for the sandbox-writable cache.
The new L13 test failed before the production change with the exact CNN error:
`cannot Fast CKKS zero-secret Add: component 0 is missing active row q4`.

New lightweight regressions use LogN10, 16 Q rows, and the stem.square-like
Scale `2^64/q14` (shared exactly by both operands). They cover:

- L13 compact ciphertext plus fully backed L13 plaintext; in-place and separate
  output; exact modular c0 oracle and unchanged c1 on every maintained row.
- Different ciphertext/plaintext/output Levels, down to Levels 1 and 2.
- Plaintext non-mutation and ciphertext non-mutation for separate output;
  metadata/Scale/Level preservation and absent output rows beyond the prefix.
- Plaintext high rows filled with invalid uint64 values or removed entirely:
  results still match the plaintext oracle and the shared core receives four rows.
- Nil/empty operands or outputs, metadata, unequal Scale, NTT/Montgomery/batching
  mismatch, missing prefix backing, nonzero c1, noncanonical residues and wrong N;
  rejection without panic or prefix-data mutation.
- Native complete-Q nonzero-c1 Add fallback with ordinary encryption/decryption;
  no Fast Add core dispatch. This fork-local test is not a genuine Standard baseline.

Maximum decoded complex error for the new compact Add Level/alias cases:
`1.10178475297e-9`, below the unchanged `1e-6` test tolerance. The tests exercise
ordinary DecryptNew/Decode, exact row-wise modular arithmetic and recorded shared
core dispatch, rather than only testing absence of an error.

Passed:

- 46 top-level tests from
  `go test ./schemes/ckks -run '^(TestPublicFast|TestFastCKKSZeroSecret|TestFastQP|TestGenericFastCKKSKeyGenerator|TestCKKSPublicKeyEncryptionWithoutFastCapability)' -count=1 -v`.
- The subsequently added native-fallback test and all new plaintext tests:
  `go test ./schemes/ckks -run '^TestPublicFastPlaintextAdd' -count=1 -v`.
- `go test ./core/rlwe ./schemes/bgv -short -count=1`.
- `go vet ./schemes/ckks ./core/rlwe ./schemes/bgv` and `git diff --check`.
- Isolated, unmodified genuine Standard source exported from
  `5dbffbdea05394de2ca3a432ed5318aa832e3f40`:
  `GOWORK=off go test ./schemes/ckks ./core/rlwe -short -count=1`.

## Stale scalar assertion and self-review

The previous compact scalar Add test expected a missing-q4 error, contradicting
the existing public scalar dispatcher and shared scalar kernel. It now checks
`2*values + 3`, decoded accuracy at 1e-6, preserved Scale/Level and compact output.
Its generic-consumer rejection value is retained using direct vector Add, which
remains unsupported for compact input and must leave the output unchanged.
The original nonzero-c1 native fallback and invalid-residue assertions remain.
No production scalar behavior or numerical threshold was changed.

Coding self-review checked all acceptance criteria, the complete diff, degree-zero
adaptation, strict prefix-only consumption, aliases, malformed input handling,
and preserved native fallback. The extra native-fallback regression passed;
no numerical defect, scope expansion or policy change was needed.

No Primary, CNN512, model, Golden, Experimental Trust document, Bootstrap or
benchmark was modified/executed. The existing measurement platform was not
extended or duplicated; validation uses maintained library test helpers.
behemoth1 can update this branch and recompile the unchanged CNN512 source.
Full CNN512/LogN16 Golden acceptance remains a remote integration task.
Status: `READY_FOR_WEB_REVIEW`.
