# Fast-CKKS specification

This document is the single detailed source of truth for Fast-CKKS work on the `fast-ckks` branch. It records the current implementation, the permanent architectural boundaries, and the future experiments. Statements marked **Target design**, **Requires source audit**, or **Future experiment** are not claims that the current code already implements them.

## 1. Project objective and priority

The application keeps ordinary CKKS API usage while the backend supplies Normal or Fast behavior:

```text
Application
    |
    | unchanged CKKS API usage
    v
Lattigo
    |
    +-- upstream/normal behavior
    |
    +-- Fast-CKKS behavior on fast-ckks branch
```

The application must not contain Normal/Fast branches. Fast behavior belongs in the backend library.

The primary engineering goal of the current stage is speed. Fast-CKKS must make CKKS execution sufficiently fast for large workloads, including thousands of CNN inferences used to collect FHE numerical-error data. If Fast execution is not substantially faster than Normal Lattigo, the current project goal has not been met. The project is intentionally insecure; security is not a current optimization constraint.

Development priority:

1. Fast execution throughput.
2. Functional correctness of the simplified Fast execution.
3. Clean extension points for later experiments.
4. Noise fidelity evaluation and improvement later.

Current Fast execution is not required to reproduce the full Standard FHE noise distribution. Do not add expensive computation merely to make Fast noise look more realistic.

## 2. Non-goals and scope

Fast-CKKS is not intended to redesign CKKS or the application, modify CNN/model code, reduce the frontend-configured parameter chain, refactor all Lattigo schemes, create Fast BFV/BGV/TFHE, or optimize unrelated code.

Primary scope:

- `schemes/ckks` and its Fast backend;
- CKKS circuits required by the application;
- `circuits/ckks/bootstrapping`; and
- directly required `core/rlwe` execution paths.

Other schemes remain untouched unless a genuinely shared primitive is unavoidably required by a CKKS call path. Preserve the original Normal implementation; Fast code remains separate where practical so both paths support correctness comparison, performance comparison, and later noise calibration.

## 3. Core design principle

> Preserve structure; elide computation.

Preserve public APIs where practical, CKKS parameter objects, metadata, structural compatibility, Normal behavior, and enough representation structure for future experiments. In Fast hot paths, do not compute expensive unneeded residue values, perform full-RNS work merely because storage exists, perform security-only work merely because Standard Lattigo does, synchronize every stored limb per operation, or invoke Standard QP/key-switch paths as a hidden fallback.

Unused Fast-state fields or residue arrays may exist structurally while their values are unmaintained, stale, or dormant. The optimization target is to avoid computing, updating, transforming, or synchronizing RNS residue values that are not needed for the Fast numerical result.

## 4. Fast constitution: logical CKKS semantics versus Fast storage

This section is the **authoritative target architecture** for future Fast-CKKS implementation. Current q0/q1/q2 code below is historical/current implementation evidence and must not override these rules.

### 4.1 Two different modulus systems

Fast must distinguish two concepts that ordinary RNS CKKS often stores in the same arrays.

**Logical CKKS modulus chain**

```text
Q = (q0, q1, ..., qL)
```

The application-provided `q_i` remain authoritative for:

- logical ciphertext Level;
- CKKS Rescale divisors;
- Scale evolution;
- public parameter identity;
- logical ciphertext congruence;
- Standard/Fast semantic comparison; and
- conversion back to an ordinary Lattigo ciphertext.

**Fast storage moduli**

```text
F = (f0, f1, f2)
```

The `f_i` are backend-private NTT-friendly primes used only to preserve a bounded lifted integer with more CRT capacity than the frontend's low logical moduli. They are not CKKS rescale primes and do not replace the logical parameter chain.

Use the names precisely:

```text
q_i                 logical CKKS modulus at index i
f_i                 Fast private storage modulus
logical residue     x mod q_i
Fast residue        X mod f_i
c0,c1,c2,...        ciphertext polynomial components
Level               highest active logical q_i index
ActiveStorageWidth  number of currently authoritative Fast storage residues
```

Do not call `f_i` a larger logical `q_i`. Do not infer logical level from ActiveStorageWidth.

### 4.2 Fast state and permanent invariants

Conceptually a Fast ciphertext state is

[
\mathcal S=(\ell,\Delta,A,X,B,D,d),
]

where:

- `ell` is the logical CKKS Level;
- `Delta` is the logical CKKS Scale;
- `A` is the active Fast storage basis, a subset of `{f0,f1,f2}`;
- `X=(X_0,...,X_d)` is the vector of coefficient-wise lifted integer polynomials represented by the Fast ciphertext components;
- `B=(B_0,...,B_d)` is a vector of proven coefficient-infinity bounds with `||X_j||_infinity <= B_j`;
- `D` records coefficient/NTT and Montgomery representation state; and
- `d` is the ciphertext degree.

Every production Fast operation must preserve all applicable invariants below.

**C1 — Logical congruence**

For each coefficient at logical Level `ell`:

[
X \equiv c_{logical}\pmod{Q_\ell},
\qquad
Q_\ell=\prod_{i=0}^{\ell}q_i.
]

Fast may use a different integer lift than Standard internally; congruence modulo the current logical modulus is the semantic requirement.

**C2 — Unique Fast reconstruction**

For active storage product

[
S_A=\prod_{f\in A}f,
]

every authoritative component coefficient must satisfy

[
\|X_j\|_\infty<\frac{S_A}{2},
\qquad 0\le j\le d.
]

Equivalently, for every proven component bound:

[
2B_j<S_A.
]

This strict centered-uniqueness condition is the correctness boundary. An engineering safety margin may be imposed in addition, but no operation may rely on an unproven wraparound.

**C3 — Residue consistency**

Every active Fast row must equal the same lift reduced modulo its storage modulus:

[
R_f = X\bmod f.
]

Dormant rows are not authoritative and must not be read as if they were synchronized.

**C4 — Logical Scale sovereignty**

Scale changes are determined only by CKKS mathematics and logical `q_i`. Fast storage primes never directly determine public Scale.

**C5 — Logical level and storage width are independent**

A logical transition

[
\ell\rightarrow\ell-1
]

does not imply

[
3\rightarrow2
]

or any other storage-width change. Storage contraction or expansion is a separate representation operation justified only by a capacity proof.

### 4.3 Fixed shared Fast storage prime set

The target software Fast backend uses one fixed three-prime storage set across all supported LogN profiles:

[
F=(f_0,f_1,f_2),
\qquad
\operatorname{bitlen}(f_i)\approx60.
]

This project prioritizes rapid software validation over per-LogN hardware specialization. Therefore the same `f_i` should be reused across supported LogN values.

Prime-selection requirements:

1. each `f_i` is an odd prime;
2. `f_i != f_j` for `i != j`;
3. each prime is NTT-friendly for the largest supported Standard-ring degree;
4. for the current maximum `LogN=16`, require at least
   [
   f_i\equiv1\pmod{2^{17}};
   ]
5. prefer `f_i < 2^60` so existing uint64 Montgomery/lazy-NTT machinery retains comfortable headroom;
6. if a path materializes an `f_i` together with logical `q_j` in one CRT basis, the combined moduli must be pairwise coprime; and
7. candidate primes must be validated through Lattigo's own prime/SubRing/NTT-constant construction before adoption.

The exact three numerical primes are an implementation choice. The constitution fixes the role and selection constraints, not particular constants.

An illustrative frontend chain

```text
logical bit widths: 53, 38, 38, 60, 40, 45
```

may therefore have the Fast compatibility view

```text
storage bit widths: 60, 60, 60, 60, 40, 45
```

without changing the logical CKKS chain. The first three entries denote widened Fast storage roles, not replacement logical rescale moduli. Rows `q3...` may still remain dormant in sparse Fast execution.

### 4.4 Active storage width

Fast has at most three authoritative private storage residues for the bounded-integer shortcut:

[
A_3=(f_0,f_1,f_2),\quad
A_2=(f_0,f_1),\quad
A_1=(f_0).
]

The backend may keep three residues even at logical Level 1 or Level 0 if the bound requires them. Conversely, it may contract when the target basis still uniquely represents the current lift.

A contraction

[
A_3\to A_2
]

is legal only if

[
|X|<\frac{f_0f_1}{2}.
]

Likewise

[
A_2\to A_1
]

requires

[
|X|<\frac{f_0}{2}.
]

Storage contraction:

- does not consume a logical Level;
- does not divide a coefficient;
- does not round;
- does not change Scale; and
- does not change plaintext semantics.

If a future operation needs more capacity, storage may be expanded again before the operation, provided the current smaller basis still uniquely reconstructs `X`; reconstruct `X` and reduce it into the added Fast modulus. Width changes are representation management, not CKKS operations.

### 4.4.1 Initial production width policy

For the first production integration of the private-F architecture, use a fixed active storage width of three:

[
A_3=(f_0,f_1,f_2).
]

This is an engineering policy, not a new correctness invariant.

The exact planner

[
w_{req}(B)=\min\{w:2B_j<S_w\ \forall j\}
]

remains authoritative for proving whether a state fits the available private-F capacity, but it does not force the runtime to contract to the minimum admissible width.

For the initial production path:

- import into width 3;
- keep width 3 across arithmetic, Rescale, ModUp, and later integration boundaries;
- do not implement automatic `3->2` or `2->1` contraction;
- do not add contraction/expansion hysteresis or adaptive-width policy yet.

A wider active basis than mathematically required is valid. Dynamic width is a later performance experiment and should be justified by wall-clock measurements against the fixed-width-3 baseline. If fixed width 3 is already sufficiently fast, contraction may remain unimplemented.

Existing width-1/2 support, required-width proofs, and exact expansion remain useful as correctness infrastructure and diagnostic capability.

### 4.4.2 Explicit stage-local width experiments

The initial production policy remains fixed width 3. This section does not enable automatic/adaptive runtime width selection.

A bounded experiment may explicitly contract an already valid private-F ciphertext from width (w) to a smaller width (w'<w) when:

[
2B_j<S_{w'}qquad\forall j.
]

Because every active row already stores the same authoritative integer (X_j) modulo its corresponding private prime, contraction to a prefix basis does not require integer reconstruction. The target rows are simply the existing prefix residues:

[
(X\bmod f_0,ldots,X\bmod f_{w'-1}).
]

The centered-capacity proof makes that prefix a unique representation of the same authoritative integer lift.

Therefore an explicit contraction may:
- allocate/copy only the prefix rows, or use an equivalent ownership-safe row truncation;
- preserve logical Level, Scale, degree, domain, metadata, and component bounds;
- avoid INTT/CRT/NTT reconstruction.

The reverse direction is different: expanding from a narrower basis to a wider basis requires reconstructing the authoritative integer (or an equivalent proven method) to generate the missing residues. For this reason, stage-local width experiments should avoid repeated (3\to2\to3\) oscillation.

Current authorized experiment:
- C2S/DFT feasibility may compare width 2 against width 3 and the existing compact LogicalQ Fast path;
- production Bootstrap remains width-3/private-F only where already accepted and otherwise keeps the current logical Fast path;
- no benchmark result automatically changes production policy. A later Primary decision is required.

### 4.5 Operation state transitions

The following table defines the target semantics.

| Operation | Logical Level | Scale | Lifted integer / bound | Fast storage |
|---|---|---|---|---|
| Add/Sub | unchanged | unchanged | `Z=X +/- Y`, conservatively `Bz <= Bx + By` | unchanged unless capacity planning requires expansion |
| Mul | unchanged | `DeltaX * DeltaY` | exact negacyclic polynomial product; conservative `Bz <= N*Bx*By` | must have enough pre-product capacity |
| Relinearize / KeySwitch | unchanged | unchanged | logical value preserved plus bounded key-switch error | unchanged unless bound requires expansion |
| Automorphism | unchanged | unchanged | permutation/sign change; coefficient infinity bound preserved | unchanged |
| Rotate | unchanged | unchanged | automorphism followed by bounded key-switch effect | unchanged unless bound requires expansion |
| Rescale | `ell -> ell-1` | `Delta/q_ell` | `Y=Round(X/q_ell)` | width change is optional and separate |
| Storage contraction/expansion | unchanged | unchanged | exactly the same `X` | representation only |
| NTT/INTT/Montgomery form | unchanged | unchanged | same mathematical polynomial | same authoritative storage moduli |
| ModUp | logical basis grows | basis-raise step itself does not redefine Scale | canonicalize the logical representative before extension | expand/repopulate Fast basis from canonical lift |
| Public/Standard boundary | unchanged | unchanged | reduce the authoritative lift modulo required logical `q_i` | Fast -> logical Q representation |

### 4.6 Add/Sub and multiplication bounds

Fast bound tracking is **per ciphertext component**, not one scalar guessed for the entire ciphertext.

For ciphertext degree `d`, maintain

[
\mathbf B=(B_0,\ldots,B_d),
\qquad
\|X_j\|_\infty\le B_j.
]

Bounds are semantic coefficient-domain facts. NTT coordinates may be numerically large; NTT/INTT does not change `B_j` because the bound belongs to the represented polynomial, not to its transform coordinates.

#### Add/Sub

For

[
Z_j=X_j\pm Y_j
]

with a missing higher-degree component interpreted as zero,

[
B_{Z,j}\le B_{X,j}+B_{Y,j}.
]

For the initial private-storage arithmetic contract, Add/Sub requires equal CKKS Scale. The output logical Level is

[
\ell_Z=\min(\ell_X,\ell_Y),
]

which is a logical restriction analogous to dropping unavailable higher logical residues; it does not change the private Fast lift. Output degree is `max(d_X,d_Y)`.

#### Negacyclic polynomial product

For two coefficient polynomials `A,B` in

[
R=\mathbb Z[X]/(X^N+1),
]

[
\|A\star B\|_\infty
\le
N\|A\|_\infty\|B\|_\infty.
]

For ciphertext multiplication,

[
Z_k
=
\sum_{i=\max(0,k-d_Y)}^{\min(d_X,k)}
X_i\star Y_{k-i}.
]

Therefore the generic safe per-component bound is

[
\boxed{
B_{Z,k}
\le
N
\sum_i
B_{X,i}B_{Y,k-i}
}.
]

For degree-one times degree-one multiplication this becomes

[
B_{Z,0}\le N B_{X,0}B_{Y,0},
]

[
B_{Z,1}\le N(B_{X,0}B_{Y,1}+B_{X,1}B_{Y,0}),
]

[
B_{Z,2}\le N B_{X,1}B_{Y,1}.
]

If only a single worst-case input bound is known, the coarser fallback is

[
\max_k B_{Z,k}
\le
N(\min(d_X,d_Y)+1)B_XB_Y.
]

This corrects the single-polynomial shortcut `N*Bx*By`, which is insufficient for ciphertext components that sum multiple polynomial products.

Raw multiplication keeps

[
\ell_Z=\min(\ell_X,\ell_Y),
\qquad
\Delta_Z=\Delta_X\Delta_Y,
\qquad
d_Z=d_X+d_Y.
]

Relinearization is a separate operation and must not be silently folded into this bound unless its own bounded error contribution has been accounted for.

#### Integer scalar multiplication

For exact integer `m`,

[
Z_j=mX_j,
\qquad
B_{Z,j}\le |m|B_{X,j}.
]

The `MulInteger` arithmetic operation itself does not consume a logical Level and does not intrinsically change Scale. If a caller deliberately multiplies both coefficients and Scale by the same integer for message-preserving normalization, the Scale update is a separate explicit metadata transition.

Measured LogN13/P93 evidence motivating the wider basis remains:

[
B_{max,current}\approx2^{121}
]

at the worst observed pre-Rescale generated-power checkpoint. The current q0q1q2 centered capacity was only about `2^133`, leaving roughly 12 bits. A first-order future estimate in which both multiplicative operands gain about 17 effective scale bits gives

[
B_{max,future}\approx2^{155}.
]

Three near-60-bit Fast storage primes provide centered capacity near `2^179`, about 24 bits above that estimate. These numbers are design evidence, not universal correctness bounds; the explicit per-operation bound gate remains authoritative.

### 4.7 KeySwitch, Relinearize and Rotate

Fast key operations must preserve the same logical Level and Scale unless the corresponding CKKS operation explicitly says otherwise.

Conceptually:

[
X' = X + E_{KS},
]

so a safe bound is

[
B'\le B+B_{KS}.
]

Current zero-secret/error-only key material is one implementation mode that makes `E_KS` bounded and cheap. The constitution does not permanently require zero secret. Any future key mode must provide a bound and Fast-storage representation that preserves C1-C4.

A pure automorphism only permutes/sign-changes coefficients and therefore preserves the coefficient infinity bound. Rotation inherits the key-switch bound after the automorphism.

### 4.8 Rescale theorem with independent storage primes

For logical Level `ell`, let

[
Q_\ell=q_\ell Q_{\ell-1}
]

and let the Fast authoritative lift satisfy

[
X=c+kQ_\ell.
]

Fast Rescale must compute

[
Y=\operatorname{Round}(X/q_\ell),
\qquad
\Delta'=\Delta/q_\ell,
\qquad
\ell'=\ell-1.
]

Because `kQ_{ell-1}` is integral,

[
\operatorname{Round}(X/q_\ell)
=
\operatorname{Round}(c/q_\ell)+kQ_{\ell-1},
]

hence

[
Y\equiv\operatorname{Round}(c/q_\ell)
\pmod{Q_{\ell-1}}.
]

Therefore the physical rounded division and the Scale update must use the **logical divisor `q_ell` even when no Fast storage row is modulo `q_ell`**.

Using `f_i` as the divisor would be a semantic error. For example, dividing coefficients by `q_ell` but Scale by `f_i` introduces a message multiplier `f_i/q_ell`.

For a proven input component bound

[
\|X_j\|_\infty\le B_j,
]

with odd logical prime `q_ell`, nearest-integer division gives the safe post-Rescale bound

[
\boxed{
B'_j
=
\left\lfloor
\frac{B_j+(q_\ell-1)/2}{q_\ell}
\right\rfloor
}.
]

This follows because ties at exactly one half cannot occur for integer numerators divided by odd `q_ell`.

After logical Rescale, optionally contract storage only when the post-Rescale bound vector proves that the smaller storage product is sufficient. Logical Rescale and storage contraction must remain separately testable operations.

### 4.9 ModUp canonicalization boundary

ModUp is different from same-level arithmetic because congruence modulo a small logical basis does not identify a unique representative in a larger logical basis.

Before raising a logical ciphertext from `Q_ell` to a larger modulus basis, Fast must first choose the canonical centered logical representative

[
C=\operatorname{Center}_{Q_\ell}(X\bmod Q_\ell),
\qquad
C\in(-Q_\ell/2,Q_\ell/2].
]

For the current Bootstrap Level-0 ModUp this reduces to

[
C=\operatorname{Center}_{q_0}(X\bmod q_0).
]

For odd \(q_0=2m+1\), the centered convention is explicit:

[
\operatorname{Center}_{q_0}(r)=
\begin{cases}
r, & 0\le r\le m,\\
r-q_0, & m+1\le r<q_0.
\end{cases}
]

Therefore the residue \(r=\lfloor q_0/2\rfloor=q_0>>1\) is the **positive** representative \(+m\). A historical implementation that branches on \(r\ge q_0>>1\) selects \(-(m+1)\) at this single residue; that value is congruent modulo \(q_0\) but is not the canonical centered representative and must not override this architecture rule.

Only then may the backend populate Fast storage residues or any newly materialized logical residues from `C`.

ModUp is therefore an explicit **logical canonicalization boundary**. Merely extending an arbitrary lift `X=c+kQ_ell` would generally produce the wrong representative in the enlarged basis.

### 4.10 Public boundary

When a Fast value must become an ordinary logical-Q Lattigo value, reconstruct or otherwise derive the authoritative lift `X` and populate each required logical row as

[
X\bmod q_i.
]

The public result must expose the frontend's original logical Level, Scale, parameter identity, and key semantics. Fast storage primes are backend-private and must not leak into public CKKS parameters.

### 4.10.1 Fusion of adjacent representation boundaries

The private-F architecture defines semantics, not a requirement to materialize every intermediate representation when no operation can observe it.

An adjacent conversion chain may be fused away when all of the following are proven:

- the input and output observable semantics are identical to the unfused private-F reference path;
- no arithmetic operation, bound transition, key operation, or caller observes the intermediate private-F state;
- the canonical lifted integer used by the fused kernel is the same authoritative integer that the unfused path would store;
- all required capacity conditions are already proven;
- no dormant/full-RNS logical rows are materialized as a consequence;
- the accepted canonicalization convention and Logical-Q metadata semantics are unchanged.

For the production Level-0 ModUp bridge, the accepted reference chain

[
q_0\text{-LogicalQ}
\to
ImportLevel0(F_3)
\to
FastStorageModUpLevel0
\to
compact\ LogicalQ
]

has no observable private-F arithmetic between import and export. Since `ImportLevel0` already selects

[
C=\operatorname{Center}_{q_0}(r)
]

and `FastStorageModUpLevel0` is idempotent on that canonical lift, a production optimization may directly compute maintained logical rows from the same `C` without physically materializing transient `F_3` residues.

This is a **fused private-F boundary kernel**, not a Standard/full-RNS fallback. The unfused private-F chain remains the semantic oracle. When later production stages actually execute arithmetic in private-F storage, those observable private-F states must remain physically distinct as required by Section 4.11.

### 4.10.2 Normalized Trace over bounded private-F lifts

For the Standard-ring Trace operation, let `g` be the exact normalization gap used by Lattigo:

[
g=2^{\operatorname{LogN}-\operatorname{logN}-1},
]

with the existing extra factor of two when `logN = 0`.

Let

[
\mathcal T_g(X)=\sum_{\sigma\in H_g}\sigma(X)
]

denote the unnormalized automorphism sum generated by the same Trace schedule.

Lattigo's Trace theorem states that a monomial either vanishes under this sum or is multiplied by `g`. By linearity, every coefficient of `\mathcal T_g(X)` is divisible by `g`. Therefore the normalized Trace has an exact integer-lift interpretation:

[
Y=\frac{\mathcal T_g(X)}{g}.
]

This is stronger than merely multiplying by `g^{-1} mod Q`: for a bounded integer lift the normalized result itself is an integer polynomial.

If

[
\|X_j\|_\infty\le B_j,
]

then each automorphism is a signed coefficient permutation, hence

[
\|\mathcal T_g(X_j)\|_\infty\le gB_j,
]

and exact coefficient-wise division gives

[
\|Y_j\|_\infty\le B_j.
]

For private-F execution, the required intermediate-capacity condition is therefore

[
\boxed{2gB_j<S_A\quad\forall j.}
]

When this holds, the entire unnormalized Trace is uniquely represented throughout the automorphism-sum tree.

In NTT private-F storage, the implementation may:

1. form the unnormalized automorphism sum directly in the F residues;
2. only after the complete sum, multiply each active F row by `g^{-1} mod f_i`.

Because divisibility by `g` is proven in the authoritative integer domain and `g` is a power of two while every `f_i` is odd, this residue-wise final normalization equals exact integer division.

Do **not** pre-multiply a private-F ciphertext by `g^{-1} mod f_i` before the automorphism sum. That modular representative need not correspond to a bounded authoritative integer lift and can violate the private-storage uniqueness invariant even though the final logical Trace would be correct modulo Q.

The normalized private-F Trace preserves:
- logical Level;
- CKKS Scale;
- degree;
- NTT domain;
- fixed production storage width 3;
- and a safe post-bound no larger than the input bound.

### 4.10.3 Private-F plaintext mirrors for CKKS linear transforms

A CKKS plaintext diagonal is not a small scalar and must not be reinterpreted from one logical residue row.

The CKKS encoder first quantizes each coefficient to a single signed integer (P_k), then stores (P_k mod q_i) in every logical CRT row. Therefore a linear-transformation plaintext has an authoritative integer-polynomial interpretation whenever that integer representative is uniquely recoverable.

For the current Fast DFT matrices, the practical one-time mirror construction is:

1. take the complete logical-Q polynomial of the encoded matrix diagonal at its declared `LevelQ`;
2. undo Montgomery form and NTT on all logical q rows;
3. CRT-reconstruct the centered coefficient
   [
   P_k=operatorname{Center}_{Q_{LevelQ}}(p_kmod Q_{LevelQ});
   ]
4. require
   [
   2|P_k|<S_3
   ]
   for every coefficient;
5. encode the same (P_k) into private storage:
   [
   (P_kmod f_0, P_kmod f_1, P_kmod f_2);
   ]
6. transform those F rows to the NTT domain for repeated execution.

This reconstruction is circuit-initialization work, not a hot-path operation. Arbitrary-precision CRT is acceptable during this one-time conversion. A future direct encoder-to-F path may replace it but must produce the same integer polynomial.

Do not:
- infer a plaintext lift from q0 alone;
- infer it from q0/q1 unless uniqueness is separately proven;
- reinterpret Montgomery/NTT q residues as F residues;
- treat matrix `Scale` as a storage modulus.

For a private ciphertext component with
[
|X_j|_inftyle B_j
]
and an integer plaintext polynomial (P), use the convolution bound
[
|X_jstar P|_infty
le
B_j|P|_1,
qquad
|P|_1=sum_k|P_k|.
]

For a diagonal linear transform
[
Y_j=sum_d operatorname{Rot}_d(X_j)star P_d,
]
automorphisms are signed coefficient permutations, so
[
oxed{
B'_{j}
le
B_jsum_d|P_d|_1.
}
]

This bound is preferred to the coarser (N B_j B_P) bound for DFT feasibility planning.

The private-F linear-transform output Scale follows ordinary CKKS semantics:
[
Delta'=DeltacdotDelta_P,
]
where (Delta_P) is the matrix plaintext Scale. Logical Level is unchanged until the explicit logical-Q Rescale transition.

The private-F implementation must prove
[
2B'_j<S_3
]
before accepting the result. If a C2S factor or factor group cannot satisfy this with the fixed width-3 production basis, it is not eligible for private-F residency without a new architecture decision.

### 4.11 Physical container separation

The target widened-storage architecture must not store residues modulo `f_i` inside coefficient rows that ordinary Lattigo code interprets as residues modulo logical `q_i`.

A row position is not merely storage: existing `ringQ.SubRings[i]`, NTT constants, Montgomery constants, Rescale code, and other logical-Q operations interpret row `i` according to `q_i`.

Therefore:

[
\boxed{\text{Fast private }f_i\text{ residues require a physically distinct storage container/basis context.}}
]

Do not create a target design in which an ordinary `rlwe.Ciphertext` row indexed as logical `q_i` silently contains data modulo `f_i`.

The target representation should keep, explicitly and separately:

```text
logical metadata:
    logical Level
    CKKS Scale
    degree / ring degree
    public parameter identity

private physical state:
    Fast storage basis F
    ActiveStorageWidth
    residue polynomials modulo f_i
    representation domain (Coeff/NTT, ordinary/Montgomery)
```

Historical q0/q1/q012 compact ciphertexts remain valid current implementation evidence because those rows are actually modulo the corresponding logical q_i. They must not be used as a template for silently relabeling widened f_i residues as q_i rows.

### 4.12 Logical-Q to Fast-storage entry boundary

Suppose an ordinary logical-Q coefficient class `c` is available at logical Level `ell`:

[
c\in\mathbb Z/Q_\ell\mathbb Z.
]

Fast import defines its boundary lift deterministically as the canonical centered representative:

[
C=\operatorname{Center}_{Q_\ell}(c)
\in(-Q_\ell/2,Q_\ell/2].
]

This is a representation choice made by the Fast boundary; it does not attempt to recover a hidden historical lift. For same-level ring arithmetic, congruence modulo `Q_ell` is the semantic requirement, and the logical Rescale theorem in Section 4.8 is invariant under adding multiples of `Q_ell`.

Import into active Fast storage product `S_A` is legal only when the chosen canonical representative fits uniquely:

[
|C|<S_A/2.
]

A generic high-level logical-Q ciphertext can therefore be rejected by the Fast importer when its canonical coefficients exceed the private storage capacity.

For Bootstrap, the preferred production boundary is **after ScaleDown reaches logical Level 0**:

[
\text{logical Q input}
\xrightarrow{\text{ScaleDown}}
q_0
\xrightarrow{\text{center/import}}
F.
]

At Level 0:

[
C=\operatorname{Center}_{q_0}(c),
\qquad
|C|\le q_0/2.
]

For the supported frontend profiles `q0` is far smaller than the three-prime Fast storage product, so this boundary has large deterministic capacity headroom. After import, Bootstrap ModUp may increase the **logical** Level while physical storage remains in the same Fast basis.

A generic import API, if retained for diagnostics or future applications, must perform the actual capacity check coefficient-wise or require an equivalent proven bound. It must never silently wrap a canonical logical representative into insufficient Fast storage.

### 4.13 Fast-storage to logical-Q exit boundary

For an authoritative unique Fast lift `X`, exporting to an ordinary logical-Q ciphertext at Level `ell` is defined coefficient-wise by:

[
r_i=X\bmod q_i,
\qquad 0\le i\le\ell.
]

This requires no equality between `X` and the canonical logical representative. Congruence is sufficient:

[
X\equiv c_{logical}\pmod{Q_\ell}.
]

The exported ordinary ciphertext must restore logical Level, Scale, degree, NTT/Montgomery state, and public parameter identity consistently.

### 4.14 Cross-basis conversion must occur through coefficient-domain integers

NTT coordinates are modulus-specific. An NTT coordinate modulo `q_i` is not the same mathematical coordinate as an NTT coordinate modulo `f_j`.

Therefore basis conversion between logical Q and private F must not copy or reinterpret NTT rows index-by-index.

Conceptually:

[
\text{source NTT residues}
\xrightarrow{\mathrm{INTT}}
\text{source coefficient residues}
\xrightarrow{\mathrm{normalize}}
\text{canonical/authoritative integer}
\xrightarrow{\bmod\ target}
\text{target coefficient residues}
\xrightarrow{\mathrm{NTT}}
\text{target NTT residues}.
]

Likewise, Montgomery-form residues must be converted to ordinary residues before integer reconstruction, and converted into the target Montgomery representation only after reduction into the target modulus.

This rule is a correctness boundary, not permission to add routine coefficient-domain round trips to hot paths. Conversion should happen only at explicit basis boundaries.

### 4.15 Storage width expansion and contraction

Storage expansion from width `w` to `w+1` is representation-only and is legal only while the smaller basis uniquely determines `X`:

[
|X|<S_w/2.
]

Then reconstruct `X` and add:

[
X\bmod f_w.
]

Contraction follows Section 4.4 and likewise changes no logical Level or Scale.

### 4.16 Fast-storage arithmetic exactness and capacity planner

This section is the frozen arithmetic contract for private-F Add/Sub/Mul foundations.

#### 4.16.1 Exactness theorem

Let

[
S_w=\prod_{i=0}^{w-1}f_i
]

for active width `w`. Suppose a Fast operation computes, in every active row, the correct modular residue of an exact integer coefficient `z`:

[
r_i=z\bmod f_i.
]

If a proven bound satisfies

[
|z|\le B,
\qquad
2B<S_w,
]

then centered CRT reconstruction from those rows returns the exact integer `z`, not merely a congruence class.

Therefore modular NTT arithmetic over the private storage primes is an exact implementation of the lifted-integer operation whenever the output bound passes the centered-capacity gate.

For ciphertexts this condition is checked for every output component bound:

[
2B_{Z,j}<S_w
\quad\text{for all }j.
]

If this proof fails, the modular rows may still be algebraically valid residues, but they no longer uniquely determine the intended lifted integer. The operation is therefore forbidden.

#### 4.16.2 Required width

Define

[
w_{req}(\mathbf B)
=
\min\left\{
w\in\{1,2,3\}:
2B_j<S_w\ \forall j
\right\}.
]

If no such width exists, the operation must fail with an explicit Fast-storage capacity error before mutating the destination.

Correctness uses the exact integer comparison above, not an approximate bit-length heuristic.

For diagnostics, implementations may additionally report a headroom metric such as

[
h
=
\log_2(S_w/2)
-
\log_2(\max(1,\max_j B_j)).
]

Headroom is observability information; the strict inequality `2B_j<S_w` is the correctness gate.

#### 4.16.3 Expand before arithmetic

For a binary operation with current widths `w_X,w_Y` and proven output bounds `B_Z`, choose an operation width satisfying

[
w_{op}
\ge
\max(w_X,w_Y,w_{req}(\mathbf B_Z)).
]

With the current monotone-width policy, use

[
\boxed{
w_{op}
=
\max(w_X,w_Y,w_{req}(\mathbf B_Z))
}.
]

Any operand narrower than `w_op` must be **exactly expanded before the arithmetic**.

Do not:

```text
compute in a too-small basis
-> wrap modulo S_w
-> expand the wrapped result
```

because once the true result exceeds the centered capacity of the narrow basis, the omitted integer lift information is lost.

Expansion is legal because the source basis already uniquely reconstructs each input lift. Conceptually:

```text
current residues
-> reconstruct authoritative X
-> X mod added f_i
-> enter the same representation domain under the widened basis
```

A future optimized RNS base-extension routine may replace explicit reconstruction only after proving equivalence.

#### 4.16.4 No implicit contraction

Add/Sub/Mul must not automatically shrink storage width.

The output width is `w_op`. Contraction remains an explicit representation-management operation governed by Sections 4.4 and 4.15.

This rule keeps arithmetic and representation policy separately testable and prevents a successful operation from hiding a subsequent narrowing decision.

#### 4.16.5 Bound provenance

A production Fast ciphertext must never carry an optimistic guessed bound.

Allowed bound sources are:

1. an exact coefficient scan at a deliberate boundary;
2. a mathematically proven operation transition from already-valid input bounds; or
3. a conservative externally proven bound encoded by the caller/specification.

At the Level-0 logical-q0 to FastStorage import boundary, every coefficient is already visited for centered lifting. The preferred initialization is therefore the exact per-component maximum:

[
B_j
=
\max_k
|\operatorname{Center}_{q_0}(c_{j,k})|.
]

This exact scan has no additional asymptotic cost at that boundary and is typically tighter than the generic `q0/2` bound.

After arithmetic, use the proven transition formulas. A later deliberate full coefficient scan may tighten a bound, but routine hot-path correctness must not depend on such scans.

Bound arithmetic itself is scalar metadata work and must be overflow-safe. Using arbitrary-precision integers for the small number of bound calculations per operation is acceptable; silent fixed-width overflow in capacity planning is forbidden.

#### 4.16.6 LogN13 capacity sanity example

For the accepted fixed storage primes, the three centered capacities are approximately:

```text
width 1: 2^59
width 2: 2^119
width 3: 2^179
```

For a 55-bit logical `q0`, a Level-0 canonical centered import has a generic component bound below about

[
2^{54}.
]

At `N=2^13`, a degree-one times degree-one multiplication has middle-component bound

[
B_{Z,1}
\le
N(B_{X,0}B_{Y,1}+B_{X,1}B_{Y,0}).
]

Under the symmetric worst-case estimate `B_{X,j},B_{Y,j} <= 2^54`:

[
B_{Z,1}
\lesssim
2\cdot2^{13}\cdot2^{108}
=
2^{122}.
]

Therefore a generic correctness proof cannot authorize width 2 for that multiplication, while width 3 has ample capacity. This is consistent with the historical measured pre-Rescale maximum near `2^121`.

This example is explanatory only. Runtime authorization still uses the exact bound vector and exact product comparison.

#### 4.16.7 Transactional capacity failure

Capacity planning and any required width expansion must complete before destructive destination writes.

If:

- `w_req > 3`;
- an operand cannot be exactly expanded;
- metadata/domain preconditions fail; or
- any bound is missing/untrusted,

the operation must return an error without leaving a partially updated Fast ciphertext that appears valid.

#### 4.16.8 Initial arithmetic-domain contract

For the first private-F arithmetic implementation:

- operands must use the same ring degree and parameter identity;
- Add/Sub require equal Scale;
- binary operations require matching coefficient/NTT domain;
- Montgomery private-F arithmetic may remain unsupported and explicitly rejected;
- raw Mul requires NTT-domain operands unless an explicit coefficient-to-NTT boundary is performed;
- logical output Level is the minimum input logical Level;
- storage width is governed only by the capacity planner above.

These restrictions are intentionally conservative. They may be relaxed later only with an explicit mathematical contract.

### 4.17 Current implementation versus this target



The current branch still implements historical q0/q1 and q0/q1/q2 maintained-residue shortcuts in which storage moduli are the frontend logical moduli themselves. That code is valid current evidence and a useful oracle for bounded behavior, but it is **not** the permanent storage architecture.

Future implementation derived from this constitution must introduce an explicit distinction between:

```text
LogicalQ[i]   = frontend q_i
FastStorage[i] = backend-private f_i
```

and must audit every place that currently reads `ringQ.SubRings[i].Modulus` to determine whether that use means logical CKKS semantics or physical Fast storage arithmetic.

## 5. Fast hot-path rules

In production Fast hot paths, avoid unless explicitly required at a deliberate boundary:

- CRT reconstruction;
- redistribution into unneeded `r2...rL` storage;
- full-Q NTT/INTT;
- full-Q arithmetic;
- QP arithmetic;
- GadgetProduct;
- Standard key switching or Standard relinearization;
- unnecessary basis extension or ModDown;
- unnecessary temporary allocations; and
- unnecessary coefficient-domain/NTT-domain round trips.

Reference/debug helpers such as r0/r1 reconstruction may remain available, but must not silently migrate into production hot paths. Dormant residue storage must not be read, synchronized, or treated as mathematically current merely because it exists. Follow Standard CKKS Level/Scale semantics when an operation requires a level transition.

### 5.1 Reconstruction representation contract

For coprime q0 and q1, `(r0, r1)` uniquely determines an element of `Z/(q0*q1)Z`. CRT uniqueness is not an architectural risk. The implementation-sensitive issue is normalization and interpretation.

Before ordinary coefficient reconstruction, the conceptual precondition is:

```text
coefficient domain
non-Montgomery representation
canonical residue r_i in [0, q_i)
```

- NTT coordinates require an appropriate INTT boundary before coefficient reconstruction.
- Montgomery values require conversion to normal representation.
- Lazy arithmetic may produce values in `[0, 2*q_i)`, `[0, 3*q_i)`, or another documented lazy range; reduce them to canonical residues before CRT interpretation.

Reference/debug code may reject nonconforming input. Production code may normalize efficiently when its caller contract permits. This requirement does not authorize routine NTT -> INTT -> NTT round trips: keep the Fast evaluator NTT-resident where possible and convert only when reconstruction is mathematically required.

### 5.2 ScaleDown semantics

For current Level `level`, define:

```text
Q_level = q0 * q1 * ... * q_level
Delta   = ciphertext Scale
rho     = Mod1 MessageRatio

currentMessageRatio = Q_level / Delta
```

Current Lattigo `ScaleDown` first removes unnecessary upper moduli while the source check is conceptually equivalent to:

```text
currentMessageRatio >= q_level * rho
```

which implies:

```text
(Q_level / q_level) / Delta >= rho
```

The loop uses `Resize`/DropLevel. It does not divide polynomial coefficients, change Scale, or introduce rounding; it only removes an unnecessary high-modulus representation. This is distinct from Rescale.

The three cases are:

1. **Already Level 0:** no level dropping and no `RescaleTo`; scalar adjustment may still occur.
2. **Free DropLevel reaches Level 0:** remove all unnecessary upper levels, apply scalar adjustment, and do not call `RescaleTo`.
3. **Free DropLevel stops at Level j > 0:** compute the current message ratio, multiply coefficients by an integer approximately `scaleUp = currentMessageRatio / targetMessageRatio`, update Scale by the same factor, then call `RescaleTo` to reach the target Bootstrap scale and eventually Level 0.

DropLevel is preferred first because it is cheaper and introduces no rescaling rounding error. If both coefficients and Scale are multiplied by `a`, the represented message is unchanged:

```text
x -> a*x
Delta -> a*Delta
(a*x)/(a*Delta) = x/Delta
```

The scalar adjustment aligns message-ratio/scale state with the Bootstrap target, whose final scale is approximately `q0 / MessageRatio` at Level 0.

### 5.3 Fast Rescale implementation rule

Section 4.8 is the authoritative mathematics. Production Fast Rescale must implement

```text
Y      = round(X / logical_q[level])
Scale' = Scale / logical_q[level]
Level' = Level - 1
```

using the authoritative Fast lift, regardless of which private storage moduli currently hold that lift.

For logical levels whose active storage row is still the same modulus as the logical divisor, an optimized ordinary drop-last implementation may remain valid. For the widened low Fast rows, where `f_i != q_i`, Standard `DivRoundByLastModulus...` cannot be used merely by dropping `f_i`: that would divide by the wrong modulus.

The widened-low-limb implementation must instead perform an exact fixed-width path conceptually equivalent to:

```text
active Fast residues
    -> unique centered lift X
    -> rounded divide by logical q[level]
    -> reduce quotient into the desired Fast storage basis
```

and then update logical Scale and Level with the same `q[level]`.

A storage-width contraction may be attempted after Rescale, but it is not part of Rescale semantics and must be guarded by the target-basis centered-capacity proof.

The existing q0/q1 and q0/q1/q2 fixed-width Rescale code is the current implementation/reference mechanism. Its assumption that maintained storage moduli and logical moduli are identical must not be copied into the target widened-storage design.

### 5.4 Target parameter profile

The supplied engineering profile, not a universal CKKS invariant, constructs 17 Q primes (`MaxLevel = 16`) in this order:

| Segment | Prime bit sizes |
|---|---|
| Q0 | 55 |
| SlotsToCoeffs Q | 39, 39, 39 |
| Circuit Q | 45 |
| EvalMod Q | eight 60-bit primes |
| CoeffsToSlots Q | four 56-bit primes |

For LogN=13, the generated profile normally has `q0 < 2^55`, `q1 < 2^39`,
and `q0*q1 < 2^94`. The N=65536/LogN=16 NTT-prime generator can select a
56-bit q0 while retaining a 39-bit q1. The supported fixed-width domain is
therefore `q0 < 2^56`, `q1 < 2^39`, and `q0*q1 < 2^95`. In both cases,
per-coefficient r0/r1 reconstruction does not require arbitrary-precision
arithmetic and fits in a two-`uint64`/`math/bits` candidate representation.
Future parameter sets must be validated rather than assumed to share this
bound.

Production direction:

```text
Reference oracle:
    existing big.Int r0/r1 reconstruction + direct rounded division

Production candidate:
    allocation-free fixed-width r0/r1 reconstruction and rounded division
    using 64/128-bit arithmetic
```

## 6. Current implementation status

This branch does not start from Phase 0. The status below is derived from the current source and branch history through `a03cec77`.

| Area | Status | Current evidence and boundary |
|---|---|---|
| r0/r1 reconstruction and redistribution utilities | Reference/debug only | `schemes/ckks/fast/fast.go`; useful as a correctness oracle or explicit materialization boundary, not the production hot path |
| First-two-limb partial NTT/INTT | Implemented | `partial_ntt.go`; transforms only residue storage for q0 and q1 |
| Fast key generation | Implemented | `keys.go`; current key lifecycle uses zero-secret semantics and Fast layouts |
| First-two-limb evaluator arithmetic | Implemented | NTT-resident Add/Sub, Mul/MulRelin, scalar/plaintext Mul, fused MulThenAdd, truncation-backed Relinearize, and Rescale operate on maintained r0/r1 residues; full `schemes.Evaluator` conformance remains blocked by its unsafe full-Q/QP `rlwe.EvaluatorProvider` surface |
| Fast multiplication | Implemented | `FastMulQ01Authoritative` avoids CRT, bounds, and redistribution; reference multiplication remains available |
| Degree truncation | Implemented | `FastTruncateDegree2To1` under the current zero-secret mode |
| Fast automorphism/rotation | Implemented | First-two-limb Standard-ring path; no evaluation-key use in the helper |
| Ring-degree conversion | Implemented | r0/r1-maintaining N1/N2 conversion |
| Fast LinearTransform | Implemented | direct and BSGS NTT/Montgomery first-two-limb paths without evaluation-key, QP, or level transition; accepts matrices encoded at a higher LevelQ |
| Fast DFT adapter | Implemented and wired | `circuits/ckks/dft/fast.go` executes CoeffsToSlots/SlotsToCoeffs factor groups with Fast LinearTransform, Conjugate, scalar arithmetic, and Rescale; bounded Stage-A Bootstrap owns the matrices and uses the same Standard scaling preparation |
| Fast polynomial execution | Implemented for the bounded single-polynomial CKKS path required by EvalMod | `circuits/ckks/polynomial/fast.go` evaluates Chebyshev polynomials with Fast q0/q1 arithmetic, reusable power/baby-step workspace, and the existing Paterson-Stockmeyer planner |
| Fast EvalMod / Mod1 | Implemented for Stage-A CosDiscrete, no-inverse, scaling=1 path | `circuits/ckks/mod1/fast.go` preserves Standard normalization, cosine offset, target-scale schedule, and DoubleAngle using Fast q0/q1 arithmetic; SinContinuous, CosContinuous, inverse/ArcSine, non-unit scaling, and iterative paths remain unsupported |
| Fast Bootstrap key support | Implemented | `circuits/ckks/bootstrapping/fast_keys.go` generates Fast evaluation-key material |
| Fast Rescale | Implemented | `schemes/ckks/fast/rescale.go` provides fixed-width q0/q1 arithmetic with Standard Scale/Level semantics; targeted tests pass |
| Fast ScaleDown | Implemented | `circuits/ckks/bootstrapping/fast_scaledown.go` preserves Standard cheap DropLevel, Level-0 integer alignment, and Fast RescaleTo semantics; it is the first stage of the bounded public Fast Bootstrap path |
| Fast ModUp basis raise / scale alignment | Implemented | `circuits/ckks/bootstrapping/fast_modup.go` restores MaxLevel structure, computes only maintained q0/q1, and feeds the Fast Trace boundary |
| Fast Trace | Implemented | `schemes/ckks/fast/trace.go` performs q0/q1-only inverse normalization, automorphism/add stages, and logN=0 final order-two stage with evaluator-owned scratch/cache |
| Fast Pack/Unpack and ring boundary | Implemented for bounded Stage-A path | `fast_packing.go` maintains q0 at Level 0 and q0/q1 at higher levels, supports same-ring and factor-two ring-degree conversion, and copies caller-owned maintained state |
| Fast compact internal ciphertext storage | Implemented for Fast-owned temporary state | `schemes/ckks/fast/ciphertext.go` preserves logical `Level` in the coefficient-row headers while allocating N-sized backing arrays only for maintained q0/q1 rows; Fast-owned internal paths avoid materializing dormant q2...qL rows, while residual Level-0/1 public active rows remain ordinary and fully materialized |
| Fast ModUp boundary | Implemented for Stage-A no-key Standard-ring path | `FastEvaluator.ModUp` performs basis raise, Fast Trace, and q0/q1-only Montgomery handoff; dense/sparse switching and key switching remain intentionally out of scope |
| Fast Bootstrap orchestration | Implemented for bounded Stage-A Standard-ring path | `FastEvaluator.Bootstrap` and `BootstrapMany` run ScaleDown → ModUp/Trace → C2S → EvalMod → S2C → unpack → public Montgomery exit without Standard/QP/key-switch fallback |

The current implementation status is not a permanent architecture claim. In particular, current zero-secret behavior and current evaluation-key material are implementation modes that future experiments may extend.

## 7. Current bounded Fast Bootstrap contract

The Stage-A public Bootstrap path currently supports only:

```text
ring.Standard
CircuitOrder = ModUpThenEncode
non-iterative parameters
EphemeralSecretWeight = 0
Mod1Type = CosDiscrete
no Mod1 inverse polynomial
degree-one NTT input in ordinary (non-Montgomery) representation
N2 = N1 or N2 = 2*N1
input Level >= 0
Residual MaxLevel <= 1
```

Fast preserves the full configured modulus chain in parameter objects, but the
active execution state maintains only q0 at Level 0 and q0/q1 at higher
levels. The public finalization converts those maintained limbs out of
Montgomery form, restores the residual default scale and metadata, and leaves
dormant higher limbs untouched. Generic Standard-Bootstrap feature parity,
iterative correction, inverse/continuous Mod1 variants, ephemeral switching,
and arbitrary ring-degree ratios are not implied.

For the bounded LogN13 full-slot profile, Fast CoeffsToSlots applies
group-specific plaintext-encoding compression to preserve exact centered
q0/q1 representability: group 0 uses exponent 4, group 1 exponent 2, and
later groups exponent 0. After each ordinary Fast Rescale, the same power of
two is restored in both maintained coefficients and metadata Scale without an
additional level. Mathematical DFT scaling and Standard DFT behavior are
unchanged; this is a bounded capacity workaround.

## 8. Current implementation versus permanent architecture

### CURRENT IMPLEMENTATION

The current Fast key path forces the secret-key contribution to zero when generating Fast material. The Fast public-key target is `(e_pk, 0)`: the error component is retained and the second component is the zero polynomial. Current evaluation-key generation is approximately error-only / `(e, 0)`-style according to the implemented Fast key-generation path. Current degree-2 truncation relies on this zero-secret mode, and current Fast operations support a bounded subset of CKKS execution.

These are current implementation facts, not permanent scientific assumptions. Do not state that `s = 0` is permanently required or that all future key material must remain zero/error-only.

### PERMANENT FAST ARCHITECTURE

The durable boundaries are: preserve Normal behavior and public structure where practical; preserve the frontend logical CKKS modulus chain as the sole authority for Level/Scale/Rescale semantics; use the fixed cross-LogN Fast private storage basis defined in Section 4 for bounded lifted-integer execution; keep logical Level independent from Fast active storage width; permit normal Level/Scale progression through Level 0; require centered-uniqueness proofs before relying on a storage basis or contracting it; canonicalize at logical ModUp boundaries; avoid hidden Standard/full-RNS fallback on stale residue storage; and leave localized extension points for alternate secrets, key models, and equivalent-noise experiments.

## 9. Key policy and future experimental knobs

Fast mode intentionally removes or elides computational contributions of relevant key material in its current implementation. The exact semantics of future nonzero-secret evaluation keys require a dedicated source-level and mathematical audit.

Planned conceptual backend experiment knobs (not implemented here) include:

```text
SecretMode:
    Zero
    StandardSampled

PublicKeyErrorDistribution
PublicKeyErrorMagnitude
RandomSeed
EquivalentNoiseModel
```

`StandardSampled` should use a normal Lattigo secret polynomial drawn from the configured `Xs`, allowing paired Normal/Fast experiments to use the same sampled secret. Do not add application/frontend Fast flags; these knobs belong in a backend experiment harness or configuration layer.

Initially, later experiments should vary Public Key error only. Do not make evaluation-key error an experiment variable yet.

For future nonzero-secret experiments, do not assume that setting the `a` polynomial to zero means all secret-dependent terms may be removed, auxiliary P moduli may be removed, or Gadget/RNS decomposition has no relevance. The roles of P, secret-dependent gadget terms, key switching, relinearization, and Galois keys require explicit audit.

## 10. Two development stages

### Stage A — Fast Engine

Current priority: complete a runnable CKKS/Bootstrap execution path with maximum practical speed. Cleaner noise than Standard FHE is acceptable. Stage A should make large workloads such as thousands of CNN inferences practical. First eliminate asymptotically or structurally expensive Standard work, then complete Fast coverage; do not prematurely micro-optimize before the end-to-end path exists.

### Stage B — Noise Experiments

After the Fast Engine works end-to-end, use the experimental knobs and, if needed, theory-guided equivalent-noise injection to study FHE numerical error. Stage B extends Stage A rather than replacing it.

Removing security-related work may alter numerical-noise propagation in key switching, relinearization, rotation, ModDown, rescaling, and bootstrapping. A later operation-specific equivalent-noise model may inject the omitted effect efficiently. Do not add guessed Gaussian noise or implement this model during Stage A.

## 11. Operations and Bootstrap audit

Before implementing additional behavior, independently trace `Mul`, `MulRelin`, `Relinearize`, `Rotate`, `Rescale`, and `Bootstrap`. Record the CKKS entry point, downstream functions, ciphertext degree, components read/written, keys used, key-switch location, RNS/level behavior, implementation layer, and smallest modification point.

Bootstrap is the highest-risk integration stage. Trace the actual current Lattigo flow before implementation, including rotations, key switching, relinearization, multiplication, rescaling, modulus/level transitions, evaluation keys, and degree changes. Do not copy the complete Standard Bootstrap implementation merely to create a Fast version; prefer the smallest reusable operation boundary.

## 12. Validation terminology and methods

Use three domains precisely:

- **Message domain:** the original vector before CKKS encoding.
- **Plaintext domain:** the encoded CKKS polynomial/plaintext representation.
- **Ciphertext domain:** RLWE ciphertext polynomial components.

Do not call the pre-encoding vector a plaintext.

### Fast correctness validation

Use Normal Lattigo as a reference where mathematically appropriate. For individual Fast primitives, compare actively maintained residues with the corresponding Normal/reference operation when contracts are expected to match. For full Bootstrap, compare decoded semantic outputs rather than requiring raw ciphertext equality unless exact equality is justified.

### Noise experiment

Given `ct_before` and `ct_after = Bootstrap(ct_before)`, define an aligned difference:

```text
ct_error = ct_after - ct_before
```

Before subtraction, ensure compatible secret, ring degree, ciphertext degree, scale, level, and domain/NTT representation. Then use:

```text
message_error = Decode(Decrypt(ct_error))
```

for statistical analysis and later CNN robustness training. Retain a sanity check against:

```text
Decode(Decrypt(ct_after)) - Decode(Decrypt(ct_before))
```

within expected CKKS tolerance.

The target workload uses polynomial ring degree `N = 65536`; this is not 65,536 CKKS slots. Standard complex CKKS normally has at most `N/2` slots. Derive slot counts from the active parameter configuration rather than hard-coding them.

## 13. Performance and testing

Performance is a first-class acceptance criterion during development. Each major Fast primitive should eventually be measurable against its Normal equivalent. Do not defer all performance measurement until semantics are complete, but do not prematurely micro-optimize before the complete Fast execution path exists. After coverage is complete, profile allocation, domain-conversion, and scratch-buffer hotspots.

Testing should include:

1. Normal regression where Normal behavior remains available.
2. Fast invariant tests for actively maintained residues, Level/Scale semantics, and intentional zero-key states.
3. Independent tests for each modified primitive.
4. Bootstrap integration and continued supported use of its result.
5. Application compatibility without Fast-specific application changes.

## 14. Known Performance Questions

Non-binding questions for investigation after Fast coverage is complete:

- Can repeated temporary polynomial allocation use reusable buffers or pools?
- Can repeated first-two-residue-limb NTT/INTT conversions be reduced?
- Can Fast ciphertexts remain in NTT domain across larger operation sequences?
- Can automorphism scratch allocations be reused?
- Can LinearTransform and Bootstrap reuse scratch buffers?

These are not requests to optimize code in this documentation task.

## 15. Hard prohibitions

- [ ] Do not modify application/model architecture or add frontend Fast/Normal flags.
- [ ] Do not reduce the configured CKKS parameter chain to two levels.
- [ ] Do not modify BFV/BGV/TFHE merely for symmetry.
- [ ] Do not infer current Lattigo behavior from names alone.
- [ ] Do not confuse logical modulus primes `q_i`, Fast private storage moduli `f_i`, logical/Fast residues, ciphertext components `c0/c1/...`, logical Level, or ActiveStorageWidth.
- [ ] Do not read or synchronize unmaintained residue storage in Fast hot paths without an explicit operation requirement.
- [ ] Do not impose an artificial Level-1 floor; follow Standard CKKS Level/Scale semantics through Level 0.
- [ ] Do not call Normal key-dependent/full-RNS paths with zero keys and label that implementation complete.
- [ ] Do not silently restore security semantics or add guessed noise during Stage A.
- [ ] Do not create repository-wide refactors.
- [ ] Do not start Bootstrap implementation before tracing its actual source path.

## 16. Definition of done

The overall project eventually satisfies:

```text
Same application source
        |
        +-- standard Lattigo -> Normal behavior
        |
        +-- Fast-CKKS fork/branch -> Fast behavior
```

The application remains unchanged. Fast-specific details remain in the backend. Stage A has a complete runnable fast path with measurable throughput advantage; Stage B can then add localized, theory-guided noise experiments without replacing the engine.
