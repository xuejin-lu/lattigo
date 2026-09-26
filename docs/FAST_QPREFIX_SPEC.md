# Fast-CKKS Q-Prefix v2 Architecture

## Status

Authoritative architecture for Secondary branch:

- repository: `xuejin-lu/lattigo`
- branch: `fast-qprefix`
- branch point: `40532b4dce5c7eeae2db5b0b6f21be64801ce923`
- branch-point milestone: `fast-ckks: validate LogN13 P93 Q012 bootstrap candidate`

This branch intentionally starts **before** the later independent private-F storage architecture.

The existing `fast-ckks` branch remains the private-F comparison branch and must not be rewritten.

---

# 1. Objective

Keep ordinary CKKS parameters, Levels, Scales, and public APIs, while making Fast execution maintain only a bounded prefix of the actual logical-Q chain.

There is no independent storage basis:

[
oxed{	ext{Fast storage residues are actual logical }q_i	ext{ residues.}}
]

The initial Q-prefix cap is index 3.

For logical Level (ell), the maintained residue set is:

[
A(ell)={q_0,ldots,q_{min(ell,3)}}.
]

Equivalently, the maintained-row count is:

[
w_Q(ell)=min(ell+1,4).
]

Therefore:

| Logical Level | Maintained Fast Q prefix |
|---:|---|
| (ellge3) | q0, q1, q2, q3 |
| 2 | q0, q1, q2 |
| 1 | q0, q1 |
| 0 | q0 |

The cap of four rows is an initial engineering policy for this branch, not a universal CKKS theorem.

---

# 2. Critical distinction from the private-F architecture

The private-F branch introduced:

[
F=(f_0,f_1,f_2)
]

whose storage width was independent of logical Level.

This branch does **not** do that.

There is no:

- `f0/f1/f2`;
- Q<->F conversion;
- independent active storage width;
- private storage prime generation;
- private-F plaintext mirror.

Every maintained row is already a genuine logical residue:

[
r_i[k]=X[k]mod q_i.
]

When logical Level decreases, the maximum available maintained prefix decreases naturally with it.

---

# 3. Logical Level remains authoritative

The public CKKS parameter chain is unchanged.

A configured chain

[
Q=(q_0,ldots,q_L)
]

remains the full configured chain.

Fast does not redefine Level.

Logical Level (ell) means the ciphertext lives modulo:

[
Q_ell=prod_{i=0}^{ell}q_i.
]

The Q-prefix rule only controls which residue rows are actively maintained in Fast hot paths.

Dormant rows above the maintained prefix may remain structurally present but stale/unallocated according to the existing compact Fast policy. They are never authoritative until explicitly materialized at a proven boundary.

---

# 4. Authoritative lifted-integer interpretation

For operations that need an integer interpretation, let the maintained prefix product be:

[
S_Q(ell)=prod_{i=0}^{min(ell,3)}q_i.
]

A component coefficient has an authoritative centered lift (X) only when its proven bound (B) satisfies:

[
oxed{2B<S_Q(ell).}
]

Under this condition, the maintained Q-prefix uniquely determines:

[
Xin(-S_Q(ell)/2,S_Q(ell)/2).
]

This centered-uniqueness rule is inherited from the later private-F work, but the storage moduli are now the real q-prefix.

Bound tracking is therefore retained as a correctness tool even though the independent F basis is removed.

If an operation cannot prove its required intermediate/result bound fits the current Q-prefix capacity, that operation is not valid under this Q-prefix policy. Do not silently read dormant q rows or fall back to Standard full-RNS.

---

# 5. Capacity profile for the current LogN13 parameter family

The historical branch-point profile has approximately:

- q0: 55-56 bits;
- q1: ~39 bits;
- q2: ~39-40 bits;
- q3: ~39 bits.

Hence the approximate prefix capacities are:

[
q_0sim2^{55},
]

[
q_0q_1sim2^{94},
]

[
q_0q_1q_2sim2^{133-135},
]

[
q_0q_1q_2q_3sim2^{172-174}.
]

The four-row product fits naturally in the existing ~192-bit fixed-width design space used by historical Q012 work.

Do not hard-code these bit lengths as universal parameter rules. Validate the actual q values from the active CKKS parameters.

---

# 6. Historical evidence at the branch point

Before private-F was introduced, the accepted LogN13/P93 candidate already used a bounded widened-Q strategy:

- ordinary Fast regions primarily maintained Q01;
- selected polynomial/Paterson-Stockmeyer regions used Q012;
- selected bounded DoubleAngle windows could use q2;
- contraction back to Q01 occurred only after centered uniqueness was proven;
- q3+ remained dormant;
- the end-to-end public Bootstrap system milestone was validated.

This historical evidence proves that an independent F basis was not required for that accepted LogN13 system milestone.

Q-prefix v2 generalizes and simplifies this idea by making the maximum maintained prefix policy explicit and uniform.

---

# 7. Arithmetic bound rules

Retain the proven bound model from the later storage research.

For each ciphertext component (j):

[
|X_j|_inftyle B_j.
]

## Add/Sub

[
B'_{j}le B_{X,j}+B_{Y,j}.
]

## Negacyclic multiplication

[
|Astar B|_infty
le
N|A|_infty|B|_infty.
]

For ciphertext multiplication:

[
B'_{k}
le
Nsum_i B_{X,i}B_{Y,k-i}.
]

## Integer scalar

[
B'_jle |m|B_j.
]

## Plaintext polynomial

For integer plaintext polynomial (P):

[
|X_jstar P|_infty
le
B_j|P|_1.
]

## Diagonal LinearTransform

[
B'_j
le
B_jsum_d|P_d|_1.
]

Every hot-path operation must check the result/intermediate bound against the actual current Q-prefix product when an authoritative lifted-integer interpretation is required.

---

# 8. Rescale semantics

CKKS Rescale remains:

[
Y=operatorname{Round}(X/q_ell),
]

[
ell'=ell-1,
qquad
Delta'=Delta/q_ell.
]

The divisor is always the logical top modulus (q_ell), even when (q_ell) is above the maintained prefix.

If the current maintained Q-prefix uniquely identifies (X), Fast may reconstruct (X) from only that prefix, divide by the logical (q_ell), and reduce (Y) into the new maintained prefix:

[
A(ell')={q_0,ldots,q_{min(ell',3)}}.
]

Use the conservative result bound:

[
B'_j=
leftlfloor
rac{B_j+(q_ell-1)/2}{q_ell}
ightfloor.
]

Require:

[
2B'_j<S_Q(ell')
]

before accepting the new state.

Important transitions:

- (ell>4	oell-1): maintained q0123 prefix unchanged;
- (4	o3): maintained q0123 prefix unchanged;
- (3	o2): output naturally contracts q0123 -> q012;
- (2	o1): q012 -> q01;
- (1	o0): q01 -> q0.

A shrinking prefix is valid only under the new-level centered-capacity proof.

---

# 9. DropLevel / logical-level reduction without Rescale

Standard DropLevel changes the logical modulus but does not divide coefficients.

For reductions that do not cross the prefix cap, no maintained-row change is necessary.

When a DropLevel crosses a maintained-prefix boundary, for example:

[
3	o2,
]

two cases exist.

## Same-lift contraction

If the existing authoritative bound already proves:

[
2B<S_Q(2)=q_0q_1q_2,
]

the q012 prefix already uniquely represents the same (X), so q3 can simply become dormant.

## Canonical logical contraction

If the same integer lift does not fit the smaller prefix but the logical operation is allowed to forget the dropped modulus, choose the new canonical logical representative:

[
C=operatorname{Center}_{Q_{ell'}}(Xmod Q_{ell'}).
]

This preserves the new logical ciphertext class and automatically fits the new full logical modulus.

Such canonicalization is an explicit semantic boundary and must not be confused with a no-op row truncation.

Which behavior applies to each production DropLevel boundary must be specified and tested; do not silently change representatives.

---

# 10. ModUp semantics

For the current Bootstrap Level-0 entry:

[
C=operatorname{Center}_{q_0}(Xmod q_0).
]

For odd (q_0=2m+1):

[
0le rle mRightarrow C=r,
]

[
m+1le r<q_0Rightarrow C=r-q_0.
]

In particular:

[
r=q_0>>1=m
]

uses the positive representative (+m).

When ModUp raises logical Level to (L), populate only the maintained Q prefix:

[
q_0,ldots,q_{min(L,3)}
]

from the same canonical integer (C).

No independent storage-basis conversion is involved.

---

# 11. Plaintext and LinearTransform

Encoded CKKS plaintext matrices already contain logical q residues.

Q-prefix execution should consume only the maintained logical rows directly.

There is no private plaintext mirror and no full-Q CRT reconstruction merely to obtain a Fast plaintext representation.

For a matrix encoded at a higher LevelQ than the current ciphertext, reuse its actual q-prefix rows consistent with the existing Fast LinearTransform contract.

This is a central engineering advantage of Q-prefix storage.

---

# 12. NTT and Montgomery domains

Q-prefix execution should preserve the existing optimized Fast representation:

- NTT-resident where possible;
- Montgomery form where the existing Fast logical path benefits from it;
- domain transitions only at mathematically required reconstruction/division boundaries.

Do not introduce a second residue representation merely for storage.

The existing evaluator-owned scratch and BSGS reuse are preferred over allocation-heavy standalone paths.

---

# 13. Compact physical storage

The logical ciphertext Level may be much larger than the number of maintained rows.

For logical Level (ell), allocate N-sized backing only for:

[
w_Q(ell)=min(ell+1,4)
]

rows in Fast-owned temporary ciphertexts where the existing compact representation permits it.

Higher logical rows remain dormant/unallocated.

At public/API boundaries, preserve the structural requirements already established by the Fast branch.

---

# 14. Safety rule for dormant rows

Dormant rows are not authoritative.

Never:

- read them in arithmetic;
- call Standard full-RNS operations that assume they are current;
- silently synchronize them on every operation;
- infer an integer lift from stale dormant rows.

If a boundary genuinely needs a previously dormant logical q row, materialize it from the proven authoritative Q-prefix lift or from another explicitly specified canonical logical representative.

---

# 15. Production policy

Initial Q-prefix v2 policy:

[
oxed{w_Q(ell)=min(ell+1,4).}
]

Do not add adaptive per-operation width selection in the first implementation.

The purpose of the cap is simplicity:

- at high Level, at most four real q rows;
- as Standard Level falls below q3, maintained rows naturally fall with it;
- no second modulus family;
- no Q<->F conversion;
- no width oscillation policy.

Later profiling may justify maintaining fewer than the cap for selected operations, but that is a separate optimization decision.

---

# 16. What is retained from the private-F research

The following ideas are retained because they are mathematical/engineering improvements independent of F:

- explicit per-component coefficient bounds;
- strict centered-uniqueness checks;
- canonical ModUp midpoint rule;
- exact logical Rescale semantics;
- transactional capacity failures;
- explicit representation/domain boundaries;
- L1-based plaintext/LinearTransform bounds;
- benchmark gates;
- separation of mathematical minimum capacity from performance policy.

The following are **not** part of Q-prefix v2:

- fixed (f_0,f_1,f_2);
- `FastCiphertext` as an independent-F production container;
- private-F plaintext mirrors;
- Q<->F conversion bridges;
- private-F resident Trace/C2S.

Those remain available only on the separate `fast-ckks` comparison branch.

---

# 17. First architecture audit before implementation

Before broad production changes, produce a reproducible LogN13/P93 table for the actual Bootstrap path.

For every relevant stage/boundary record:

- stage name;
- logical Level;
- maintained prefix dictated by (w_Q(ell));
- actual q values and prefix product;
- proven/observed component bounds;
- strict capacity result;
- operation about to execute;
- whether dormant rows are required;
- whether the existing historical implementation already supports the needed prefix.

At minimum cover:

- ScaleDown;
- Level-0 ModUp;
- Trace;
- all C2S groups;
- all EvalMod polynomial/PS stages;
- DoubleAngle;
- all S2C groups;
- finalization/packing boundaries.

The primary question is:

[
oxed{
	ext{Does every required authoritative lift fit }
Q_{min(ell,3)}
	ext{ at its actual logical Level?}
}
]

If yes, an independent private-F basis is unnecessary for this profile.

If no, report the first failing stage and the exact capacity deficit before proposing any new storage basis.

---

# 18. Performance principle

The Q-prefix branch exists to test whether reusing actual logical residues is simpler and faster than maintaining a second storage basis.

Performance comparisons should use:

- current historical Q-prefix Fast baseline at the branch point;
- later private-F branch results as external comparison evidence;
- same LogN13 profile/workload where possible.

Do not sacrifice the existing highly optimized q-prefix evaluator scratch/NTT/Montgomery paths merely to mimic the private-F implementation structure.

---

# 19. Hard prohibitions

Do not:

- add independent F storage primes on this branch;
- reset or rewrite the existing `fast-ckks` branch;
- reduce the frontend CKKS parameter chain;
- change Standard Lattigo semantics;
- treat q3 as permanently active below logical Level 3;
- force four maintained rows at Level 0/1/2;
- read dormant residues as authoritative;
- hide a capacity failure with full-RNS fallback;
- weaken the accepted public Bootstrap error gates merely to make Q-prefix fit;
- start a broad code rewrite before the architecture audit.

---

# 20. Decision criterion

The independent-F architecture is unnecessary for the current target profile if the audit proves all required production states satisfy:

[
2B_j<S_Q(ell)
]

under the capped Q-prefix policy, and the resulting Q-prefix implementation preserves the required Fast correctness and performance goals.

A failure should be reported as a concrete stage/capacity counterexample, not as a general claim that Q-prefix cannot work.
