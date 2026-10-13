# Real AME(4,6) states: families, symmetry, and the classification of the real points

[![Software Heritage](https://archive.softwareheritage.org/badge/swh:1:rel:16ff6cb9d2a327de3a649e4c0fa36485de63f5b3/)](https://archive.softwareheritage.org/swh:1:rel:16ff6cb9d2a327de3a649e4c0fa36485de63f5b3)
[![release](https://img.shields.io/github/v/release/palashpandya/real-ame46)](https://github.com/palashpandya/real-ame46/releases/tag/v1.0)
[![code: MIT](https://img.shields.io/badge/code-MIT-blue)](LICENSE)
[![data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey)](LICENSE-DATA.md)

Companion code and data for the article

> **Euler's thirty-six officers need not all be entangled: real AME(4,6), families and classification**
> Palash Pandya (in preparation, 2026; the reference will be added here once it is available)

**Author.** Palash Pandya
Center for Theoretical Physics, Polish Academy of Sciences, Aleja Lotników 32/46, 02-668 Warsaw, Poland
Center for Quantum-Enabled Computing, Center for Theoretical Physics, Polish Academy of Sciences, Warsaw, Poland
Contact: ppandya@cft.edu.pl · GitHub [@palashpandya](https://github.com/palashpandya)

An AME(4,6) state is a pure state of four six-level systems whose two-party marginals are all maximally mixed;
equivalently, six times any of its 36 × 36 reshapings is a *2-unitary* matrix (the quantum version of Euler's
36 officers). Every AME(4,6) state reported before this work has complex coefficients. Rather et al.
(PRL 128, 080507, 2022, Supplemental Material) conjectured that no **real orthogonal** 2-unitary of order 36 exists,
and the 2026 review (Rep. Prog. Phys. 89, 057601, open problem T2) asks: *is there an AME(4,6) state with real
coefficients?*

This repository contains, with self-contained verification scripts:

* an explicit **real** AME(4,6) state ψ with 60 non-zero coefficients in {±1, √2}/√72 (exact proof in ℚ(√2)),
  whose three 2-unitaries are direct sums of 1 × 1 and 2 × 2 orthogonal blocks;
* the **local symmetry** behind it: ψ lies in the −1 eigenspace of u^⊗4, u = diag(i, −1, −1, −i, −i, −i), the
  sector in which it was discovered; on its support the AME conditions become linear;
* two **continuous families** of AME(4,6) states through ψ, in closed form: a 10-parameter family F₇₂^± on 72 basis
  states (two complex-conjugate components, with uniform-amplitude members that are 2-unitary complex weighing
  matrices) and the 11-parameter faces F₆₀ and F₆₀′ on 60 basis states at its two ends;
* the **classification** of the real AME states on both faces: 4096 states (up to local phases), 1024 of them limits
  of F₇₂^±, in **exactly 20 classes** under local unitaries and party permutations, 8 on F₆₀ and 12 on F₆₀′,
  with a rigorous (interval-arithmetic) proof that the 20 classes are pairwise inequivalent.

## The state

√72 |ψ⟩ = Σ_ab |ab⟩ ⊗ |φ_ab⟩, with the officer states |φ_ab⟩ (classical officers in bold):

| φ_ab | b=0 | b=1 | b=2 | b=3 | b=4 | b=5 |
|---|---|---|---|---|---|---|
| a=0 | −\|12⟩+\|21⟩ | −\|02⟩+\|20⟩ | −\|01⟩+\|10⟩ | **√2\|44⟩** | **√2\|55⟩** | **√2\|33⟩** |
| a=1 | −\|01⟩−\|10⟩ | \|00⟩+\|54⟩ | \|35⟩+\|43⟩ | −\|15⟩−\|32⟩ | −\|23⟩+\|41⟩ | −\|24⟩+\|52⟩ |
| a=2 | \|02⟩+\|20⟩ | −\|35⟩+\|43⟩ | \|00⟩−\|54⟩ | −\|25⟩+\|31⟩ | −\|13⟩+\|42⟩ | \|14⟩+\|51⟩ |
| a=3 | **√2\|34⟩** | −\|23⟩−\|41⟩ | \|13⟩+\|42⟩ | **√2\|50⟩** | −\|11⟩+\|22⟩ | **√2\|05⟩** |
| a=4 | **√2\|53⟩** | \|15⟩−\|32⟩ | \|25⟩+\|31⟩ | \|11⟩+\|22⟩ | **√2\|04⟩** | **√2\|40⟩** |
| a=5 | **√2\|45⟩** | −\|14⟩+\|51⟩ | \|24⟩+\|52⟩ | **√2\|03⟩** | **√2\|30⟩** | \|12⟩+\|21⟩ |

The 12 classical officers give the set P of basis states with coefficient 1/6; the 24 entangled officers give the
set E of 48 basis states with coefficient ±1/√72. In each of the three reshapings, O = 6M is, up to permutations of
rows and columns, a direct sum of twelve 1 × 1 blocks (1) and twelve blocks (1/√2)[[s₁, s₂], [s₃, s₄]] with
sᵢ = ±1 and s₁s₃ + s₂s₄ = 0. So O is real orthogonal and 2-unitary, and the proof (Proposition 1 of the article)
can be checked by hand. In the officers picture, each reshaping has 12 classical officers and 24 officers in
two-term Bell-type states.

## Results and how they are established

Section, proposition, table and figure numbers refer to the current draft in `paper/`.

| # | Statement | Paper | Script | Status |
|---|---|---|---|---|
| 1 | ψ is AME(4,6): all 1- and 2-party marginals are maximally mixed | Prop. 1 | `01_state_exact.py` | exact, ℚ(√2) |
| 2 | Block structure above; 12 classical + 24 Bell-type officers in every reshaping (rows and columns); 12 product states with overlap 1/6, the maximum for any AME(4,6) state | Sec. II, Fig. 1 | `01` | exact |
| 3 | The 72 states E ∪ P ∪ P′ form an orthogonal array OA(72,4,6,2) of index 2; per cut 12 E-blocks and 6 mixed blocks with P on one diagonal and P′ on the other. The 10 (resp. 11) parameters of the families are independent modulo local unitaries; at ψ the 14 deformations are 11 face directions + 3 spanned by the branch tangents | Secs. III, IV | `02_families.py` | exact; independence and 11 + 3 numerical |
| 4 | **Support-72 family.** Up to diagonal local unitaries, the AME states with support exactly E ∪ P ∪ P′ are the states Ψ_δ(α, β), \|δ\| < 1, that satisfy the mixed congruences (see [closed forms](#the-closed-forms)). They form two components F₇₂^± ∋ Ψ_δ(0, ±π/2), exchanged by complex conjugation, with 10 parameters each | Prop. 2(a) | `02`, `14_family_proof.py` | exact (Smith normal form) + numerical spot checks |
| 5 | **Faces.** The AME states with support exactly E ∪ P are the states Ψ₁(α) (face F₆₀, 11 parameters, one component, moduli rigid, ψ = Ψ₁(0)); those with support E ∪ P′ are the Ψ₋₁(β) (face F₆₀′, ψ₋ = Ψ₋₁(0)) | Prop. 2(b) | `02`, `14` | exact + numerical |
| 6 | The closures of F₇₂^± meet F₆₀ exactly in the 9-torus C on which the sums of α over the quartets P₁, P₂, P₃ agree, and F₆₀′ in the analogous C′. **Eight straight branches** δ ↦ Ψ_δ(0, β), β = η_k π/2 on P′_k (η ∈ {±1}³), end at ψ, four in each component | Prop. 2(c), Fig. 3 | `02`, `14` | exact |
| 7 | At δ = −1 each branch from ψ ends at Ψ₋₁(β) with β = ±π/2 on each quartet, a **real** state of the opposite face up to a common phase on P′; for constant β this is ψ₋ = (Σ_E ε\|x⟩)/√72 + (Σ_P′ \|x⟩)/6 | Sec. IV | `02` | exact |
| 8 | No state of F₇₂^± is real, even up to diagonal local unitaries: the invariant ω ≡ ±π/2 separates the components (and the mixed congruences have no solution over GF(2)). At δ = 0 the states have uniform amplitudes: √2·6M has entries in {0, ±1, ±i} with W W† = 2I, a 2-unitary complex weighing matrix. **Uniform states:** up to diagonal local unitaries F₇₂⁺ contains a uniform state with coefficients ±1/√72 except ±i/√72 on the quartet P′₁, and a straight path of AME states joins it to ψ; four imaginary coefficients are the minimum for coefficients in {±1, ±i}/√72 (every block has an even number, ω needs an odd number on the 18 states of P ∪ P′ with x₁ ≠ 0), and the minimal sets are the six quartets | Prop. 2(d), Sec. IV, Figs. 2, 3 | `02`, `14` | exact; path numerical |
| 9 | **Proof of Proposition 2**, every computational fact it uses: the block structure on the three supports; the 108 E-moduli equations have rank 48; the 36 E-congruences have rank 28, all invariant factors 1, and are solved by the diagonal local phases (rank 20); the quartet structure of the mixed blocks; the 18 mixed congruences have rank 14 and invariant factors 1¹³ 2; ω is gauge invariant with ω(Ψ(0, ±π/2)) = ±π/2; the E-fixing diagonal gauge is a common phase and the sign change of P₃ ∪ P′; equal quartet sums of α give exactly 8 limiting β (4 per component), unequal sums none; the uniform state with 4 imaginary coefficients and the minimality of 4 | Prop. 2 (proof) | `14` | exact + numerical spot check |
| 10 | **Symmetry:** all 72 support points lie in the 360-dim −1 eigenspace of u^⊗4, u = diag(i, −1, −1, −i, −i, −i), so every state of both families (ψ included) satisfies u^⊗4\|φ⟩ = −\|φ⟩. Diagonal symmetry group of the support: Z₂ × Z₂ × Z₄ | Sec. III | `03_symmetry.py` | exact |
| 11 | On each face, 2³² real sign patterns give **2048 real AME states** modulo local phases (512 on the gluing locus). The support automorphism group (party permutations × level relabellings, officers to officers) has order 192, with 16 orbits on F₆₀ (8 × 64, 8 × 192) and 24 on F₆₀′ (4 × 32, 20 × 96). No monomial map relates the two faces | Sec. V | `04_real_points.py` | exact |
| 12 | The real states are **isolated** modulo local orthogonal maps. Most orbits have no real first-order deformation; at the others (2 + 4 orbits), the 2-dim space of real first-order deformations is obstructed at second order | Sec. V | `04` | numerical (obstruction ≥ 0.35, noise ~10⁻¹⁵) |
| 13 | **Exactly 20 classes** under local unitaries and party permutations (8 + 12). Upper bound: H₁₂^⊗4 plus a monomial map merges orbits of the same face. Lower bound: the pair (n(149/8), n(3/2)), n(t) = number of eigenvalues below t of 72·JᵀJ (J = Jacobian of the six 2-party marginals), takes 20 distinct values | Thm. 1, Table II | `05_lu_classes.py`, `07_certify_classes.py` | upper bound exact; lower bound **rigorous** (Arb interval arithmetic + Sylvester inertia) |
| 14 | First-order deformation dimensions (LU invariant): ψ and the real classes on the gluing locus: 14; the other real classes and generic face points: 11; F₇₂^±: 10; golden state: 4; Bruzda–Życzkowski H(0): 3. So no state here is LU-equivalent to those | Secs. IV, V | `06_deformations.py` | numerical (gap 10⁻¹⁵ vs 10⁻²); for the real classes the dimension is certified in `07` (as n(1/2) = 141 + dim) |
| 15 | **Where the real states lie.** Each face is an 11-torus (origin ψ, resp. ψ₋) and its real states are exactly its 2¹¹ points of order ≤ 2. The closures of F₇₂^± meet the faces in the 9-tori C, C′, whose 2⁹ = 512 order-2 points are the real limits; the other 3072 real states lie only on the faces. Branches from real points of C end at real states only for class 1 (at classes 9 and 10); the eight branches through ψ end, up to diagonal local unitaries, at ψ₋ and at ψ₋ with P′₁, P′₂ or P′₃ negated (classes 9, 9, 10, 9), one branch from each component at each. In the code coordinates θ (below) the face tori are ℝ¹¹/Λ with Λ/2πℤ¹¹ ≅ ℤ₂³ (F₆₀; hence θⱼ ∈ {0, π} reach only 256) and ℤ₂² (F₆₀′) | Sec. V | `12_real_locations.py` | exact |
| 16 | Discovery, step 0: of the order ≤ 6 local-symmetry classes, 204 pass the spectral test and 121 also the population LP; the order-4 class of u (eigenvalue −1, dimension 360) is one of them | Sec. III | `09_sector_filters.py` | exact |
| 17 | Discovery, step 1: a random-start search in the 360-dim sector of u^⊗4 finds AME states within a few restarts; most hits have deformation dimension 10, like F₇₂^± (about 1 in 8 lie on a different, 3-parameter family) | Secs. III, VI | `08_sector_search.py` | stochastic demonstration |
| 18 | Discovery, step 2: smoothed-L1 minimisation over commutant frames brings a dimension-10 hit to exactly 72 coefficients; after a Gauss–Newton polish it is an exact AME state whose support is E ∪ P ∪ P′ up to relabelling, with the moduli of the family | Sec. III | `10_sparsify.py` | numerical (seeded) |
| 19 | The tables in `paper/inputs/` and the class table agree with the data: the state (Table I), the quartets of P and P′, the class table (Table II), and the K and branch tables of the θ coordinates; the complete block lists of the three reshapings are printed by `01` | Tables I, II | `11_paper_tables.py`, `01` | exact |
| 20 | **Entangling power.** As two-quhex gates, 6ψ and every real state (all 4096, each of the 24 orderings of the parties) have the maximal entangling power e_p = (d−1)/(d+1) = 5/7, i.e. normalised e_p = 1: E(U) = E(US) = E(S) = 35/36. So each is a real orthogonal gate of maximal entangling power. Computed from Zanardi's formula and, independently, as the exact Haar average over product inputs; on the 36 basis inputs ψ gives 12 product outputs and 24 with linear entropy 1/2 | – | `13_entangling_power.py` | numerical (exact via the orthogonality proved in `01`) |
| 21 | **Klein invariant K** of Hontarenko and Życzkowski, the degree-8 contraction of four copies of T = 6ψ with four conjugate copies whose indices are matched by the Klein four-group (invariant under local unitaries and party permutations): on the 20 class representatives K = 72 (classes 1–2), 40 (3–8), and −16, −48, −32, 0, −40, −72, −56, −56, −88, −40, −40, −56 (9–20), ten distinct values. The value 108 reported for their 96-term real state is none of them, so that state is not equivalent to any of the 4096 real states | Sec. VI | `15_klein_invariant.py` | exact, ℤ[√2]; invariance numerical |

**Real classes** (Table II; data/real_classes.json). Face F₆₀ (key "+" in the data file): representative = ψ with the listed coefficients 1/6
negated. Face F₆₀′ (key "-"): representative = ψ₋ with the listed coefficients 1/6 negated. "Limit" marks the classes that
are limits of F₇₂^± (i.e. lie on C or C′).

| class | face | flipped | states | automorphism orbits | limit | deformation dim | n(149/8) | n(3/2) |
|---|---|---|---|---|---|---|---|---|
| 1 | F₆₀ | – | 128 | 64 + 64 | yes | 14 | 1048 | 157 |
| 2 | F₆₀ | 0344, 3505 | 384 | 192 + 192 | yes | 14 | 1046 | 157 |
| 3 | F₆₀ | 5303 | 128 | 64 + 64 | no | 11 | 1050 | 153 |
| 4 | F₆₀ | 5430 | 128 | 64 + 64 | no | 11 | 1044 | 152 |
| 5 | F₆₀ | 5045 | 128 | 64 + 64 | no | 11 | 1053 | 153 |
| 6 | F₆₀ | 0533 | 384 | 192 + 192 | no | 11 | 1055 | 153 |
| 7 | F₆₀ | 0455 | 384 | 192 + 192 | no | 11 | 1052 | 153 |
| 8 | F₆₀ | 0344 | 384 | 192 + 192 | no | 11 | 1052 | 152 |
| 9 | F₆₀′ | – | 288 | 96 + 96 + 96 | yes | 14 | 1021 | 157 |
| 10 | F₆₀′ | 0434, 3055 | 96 | 32 + 32 + 32 | yes | 14 | 1029 | 173 |
| 11 | F₆₀′ | 0434, 5340 | 96 | 96 | yes | 14 | 1025 | 159 |
| 12 | F₆₀′ | 0353, 4044 | 32 | 32 | yes | 14 | 985 | 155 |
| 13 | F₆₀′ | 0353 | 288 | 96 + 96 + 96 | no | 11 | 1004 | 152 |
| 14 | F₆₀′ | 0434 | 288 | 96 + 96 + 96 | no | 11 | 1032 | 160 |
| 15 | F₆₀′ | 0545 | 288 | 96 + 96 + 96 | no | 11 | 1022 | 152 |
| 16 | F₆₀′ | 3055 | 288 | 96 + 96 + 96 | no | 11 | 1026 | 152 |
| 17 | F₆₀′ | 0434, 5405 | 96 | 96 | no | 11 | 1038 | 168 |
| 18 | F₆₀′ | 5033 | 96 | 96 | no | 11 | 1020 | 152 |
| 19 | F₆₀′ | 5340 | 96 | 96 | no | 11 | 1012 | 152 |
| 20 | F₆₀′ | 5405 | 96 | 96 | no | 11 | 1030 | 152 |

## The closed forms

**Linearisation (Lemma 1).** On the supports below, every reshaping of 6ψ is, after permuting rows and columns, a
direct sum of 1 × 1 and 2 × 2 blocks. Writing ψ_x = √(s_x/72) e^{iφ_x}, the state is AME if and only if s_x = 2 on
every 1 × 1 block and every 2 × 2 block with rows r, r′ and columns c, c′ satisfies

```
s_rc = s_r'c',   s_rc' = s_r'c,   s_rc + s_rc' = 2,     φ_rc − φ_rc' − φ_r'c + φ_r'c' ≡ π  (mod 2π).
```

The moduli and the phases decouple, and real solutions are the solutions of the same congruences over GF(2).

**The families (Proposition 2).** With ε_x = ±1 the sign of ψ_x on E, δ ∈ [−1, 1], phases α_x on P and β_x on P′:

```
Ψ_δ(α, β) = (1/√72) [ Σ_E ε_x |x>  +  √(1+δ) Σ_P e^{iα_x} |x>  +  √(1−δ) Σ_P′ e^{iβ_x} |x> ]

P₁  = {0344, 3505, 4053, 5430}    P′₁ = {0353, 3530, 4044, 5405}
P₂  = {0455, 3034, 4540, 5303}    P′₂ = {0434, 3055, 4503, 5340}
P₃  = {0533, 3350, 4404, 5045}    P′₃ = {0545, 3304, 4450, 5033}
```

On E ∪ P ∪ P′ each reshaping consists of 12 E-blocks, which carry the coefficients of ψ, and 6 *mixed* blocks with
two states a, b ∈ P on one diagonal and two states a′, b′ ∈ P′ on the other. The 18 pairs {a, b} are the pairs
within the quartets P_k, and likewise for P′. For |δ| < 1 the only remaining AME conditions are the mixed congruences

```
α_a + α_b ≡ β_a' + β_b' + π   (mod 2π),     one for each of the 18 mixed blocks,
```

whose solutions form two 10-dimensional tori, the components F₇₂^±. At δ = 1 the coefficients on P′ vanish and α is
free: the face F₆₀ = {Ψ₁(α)}, with ψ = Ψ₁(0). At δ = −1 the face F₆₀′ = {Ψ₋₁(β)}, with ψ₋ = Ψ₋₁(0).

**In the code.** Scripts 01–13, `data/family.json` and `wolfram/` use the coordinates in which the family was first
written down, with an integer 72 × 9 phase matrix K (`data/family.json`; the data files call E, P, P′ L, O, N):

```
Ψ_σ(δ, θ) = (1/√72) [ Σ_E ε e^{i(Kθ)_x} |x>  +  √(1+δ) Σ_P e^{i(Kθ)_x} |x>  +  iσ √(1−δ) Σ_P′ e^{i(Kθ)_x} |x> ],   θ ∈ ℝ⁹, σ = ±1
Φ(θ, θ₁₀, θ₁₁) = (1/√72) Σ_E ε e^{i(Kθ)_x} |x>  +  (1/6) Σ_P e^{i(Kθ)_x + iθ₁₀[x=3505] + iθ₁₁[x=5045]} |x>
```

Here K also puts phases on E; by Proposition 2 each of these states equals some Ψ_δ(α, β) up to a diagonal local
unitary, Ψ_σ(δ, θ) lies in F₇₂^σ (at θ = 0 it is Ψ_δ(0, σπ/2)), and ψ = Φ(0) = Ψ_±(1, 0). Script 14 works directly in
the paper's (α, β) form.

`wolfram/verify_families.wls` checks symbolically, in Mathematica, that both closed forms in the θ coordinates are AME
for all parameter values (δ = cos 2a).

## How the state was found

The input was only N = 4 and d = 6. We searched for AME states inside eigenspaces of local symmetries
u^⊗4 with u = V D V†, where V is Haar random and D is diagonal of finite order with a spectrum that passes two
simple necessary conditions (121 classes). One sector, with spectrum {i, −1, −1, −i, −i, −i} on every party and
eigenvalue −1 (dimension 360), produced an AME state. In the eigenbasis of the symmetry it has support 348.
Minimising the support over local frames inside the commutant U(1) × U(2) × U(3) of u (and over U(6)) gave support
72, an orthogonal array of strength 2 and index 2. On it the AME conditions are linear, which gives the family in
closed form. At δ = ±1 twelve entries vanish, and the remaining congruences have a solution over GF(2): ψ.
Script 09 reproduces the choice of the sector, 08 the search in it, and 10 the step from a search hit to the
support-72 family.

## Running

Python ≥ 3.10 with `numpy` and `scipy`; `07_certify_classes.py` also needs `python-flint` (Arb). The figure scripts in
`paper/` also need `matplotlib`; the symbolic check needs Mathematica (`wolframscript`).

```bash
pip install -r requirements.txt
./run_all.sh                          # runs 01–07 and 09–15, writes results/*.log (about 7 minutes)
python scripts/08_sector_search.py 30 0                 # optional: the stochastic search demonstration
python wolfram/make_wls.py && wolframscript -file wolfram/verify_families.wls   # optional: symbolic check
```

Scripts 01–07 and 09–15 print `[PASS]/[FAIL]` lines and exit non-zero on failure. `results/` holds the logs of a
reference run (including `08` and the Mathematica check), and `results/environment.txt` records the Python and
package versions behind them (`run_all.sh` rewrites it on every run).

## Layout

```
data/psi_terms.json               exact coefficients of √72 ψ ("+1", "-1", "+sqrt2")
data/psi.npy                      ψ as a float64 vector (index = 216a + 36b + 6c + d)
data/family.json                  the 72 points (type L/O/N = E/P/P′ of the paper, sign ε, row of K), face columns, branch table, θ′
data/uniform_min_imag_terms.json  a uniform-amplitude member with only 4 imaginary entries (i^k / √72)
data/real_classes.json            the 20 classes of real states (both faces)
data/reference/                   golden 2-unitary of Rather et al. (numerical, support 112);
                                  matrix B of Bruzda–Życzkowski, H(0) = exp(iπB/3)/36
ame46real/                        core (cuts, AME test, blocks), lattice (Smith form, GF(2)), family (closed forms),
                                  jacobian (deformations, LU-invariant spectrum), realpoints (real points, automorphisms of both faces),
                                  sector (search in the u-sector, support minimisation over commutant frames)
scripts/01..15                    the verifications listed above
wolfram/make_wls.py               generates verify_families.wls from data/family.json
wolfram/verify_families.wls       symbolic check of the closed forms
results/                          logs of a reference run
paper/                            the manuscript (main.tex, references.bib, inputs/, figures/);
                                  make_state_table.py writes inputs/state_display.tex (Table I) from data/,
                                  make_amplitudes72.py writes figures/amplitudes72.pdf (Fig. 2); scripts/11 checks the tables
```

## How to cite

Please cite the article (reference to be added on publication) and, for the code and data:

```bibtex
@misc{pandya2026realame46,
  author       = {Palash Pandya},
  title        = {Companion code and data: real {AME}(4,6) states, version 1.0},
  howpublished = {\url{https://github.com/palashpandya/real-ame46/tree/v1.0}},
  note         = {archived as \url{https://archive.softwareheritage.org/swh:1:rel:16ff6cb9d2a327de3a649e4c0fa36485de63f5b3}},
  year         = {2026}
}
```

Version 1.0 (tag `v1.0`, commit `1f591ae`) is archived in Software Heritage. Its identifiers are
`swh:1:rel:16ff6cb9d2a327de3a649e4c0fa36485de63f5b3` (the tag) and
`swh:1:dir:6695ce9630e5ced8ad497862eb57ba0680389319` (its directory). Either one, appended to
`https://archive.softwareheritage.org/`, opens the archived copy. The same identifiers are listed on the
[GitHub release](https://github.com/palashpandya/real-ame46/releases/tag/v1.0).

## Licence

Code (`ame46real/`, `scripts/`, `wolfram/`, `run_all.sh`): MIT, see `LICENSE`. Data and results: CC BY 4.0,
see `LICENSE-DATA.md` (which also lists the sources of `data/reference/`). The manuscript in `paper/` is not covered.
Copyright (c) 2026 Palash Pandya.

## AI assistance

Anthropic's Claude models (Opus 5.5 and Fable 5.1) were used as assistants throughout this work: for writing and
checking the code in this repository, for generating the tables and figures from the data, and for drafting and
editing the paper and this README. Every statement in the paper is verified by the scripts here, and the author has
reviewed all content and takes full responsibility for it.

## References

- D. Goyeneche et al., PRA 92, 032316 (2015): AME states and multi-unitary matrices.
- P. Zanardi, C. Zalka, L. Faoro, PRA 62, 030301 (2000); P. Zanardi, PRA 63, 040304 (2001): entangling power and its
  expression through operator entanglement.
- S. A. Rather, S. Aravinda, A. Lakshminarayan, PRL 125, 070501 (2020): 2-unitaries are the gates of maximal entangling power.
- S. A. Rather et al., PRL 128, 080507 (2022): the golden AME(4,6) state; conjecture in the Supplemental Material.
- G. Rajchel-Mieldzioć et al., Rep. Prog. Phys. 89, 057601 (2026), arXiv:2508.04777: review, open problem (T2).
- W. Bruzda, K. Życzkowski, arXiv:2401.01671: two-unitary complex Hadamard matrices of order 36 (H(0)).
- A. S. Hedayat, N. J. A. Sloane, J. Stufken, *Orthogonal Arrays: Theory and Applications* (Springer, 1999).
- F. Johansson, IEEE Trans. Comput. 66, 1281 (2017): Arb, the ball arithmetic used by `07_certify_classes.py`.
- M. Hontarenko, K. Życzkowski, unpublished (announced October 2026): the invariant K of `15_klein_invariant.py`.
