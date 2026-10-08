# Real AME(4,6) states: families, symmetry, and the classification of the real points

Companion code and data for the article *(title and authors to be added)*.

An AME(4,6) state is a pure state of four six-level systems whose two-party marginals are all maximally mixed;
equivalently, six times any of its 36 × 36 reshapings is a *2-unitary* matrix (the quantum version of Euler's
36 officers). Every AME(4,6) state reported so far has complex coefficients. Rather et al.
(PRL 128, 080507, 2022, Supplemental Material) conjectured that no **real orthogonal** 2-unitary of order 36 exists,
and the 2026 review (Rep. Prog. Phys. 89, 057601, open problem T2) asks: *is there an AME(4,6) state with real
coefficients?*

This repository contains, with self-contained verification scripts:

* an explicit **real** AME(4,6) state ψ with 60 non-zero coefficients in {±1, √2}/√72 (exact proof in ℚ(√2));
* two **continuous families** of AME(4,6) states through ψ, in closed form: a 10-parameter family on 72 basis
  states (with uniform-amplitude members) and an 11-parameter family on the 60 basis states of ψ;
* the **local symmetry** they all share: every state of both families lies in the −1 eigenspace of
  u^⊗4, u = diag(i, −1, −1, −i, −i, −i), the sector in which the family was discovered;
* the **classification** of all real AME states on the two support-60 faces of the family (the support of ψ and
  the opposite face): 4096 states (up to local phases) in **exactly 20 classes** under local unitaries and party
  permutations, 8 on ψ's face and 12 on the other, with a rigorous (interval-arithmetic) proof that the 20 classes
  are pairwise inequivalent.

## The state

```
√72 |ψ> = −|0012>+|0021>−|0102>+|0120>−|0201>+|0210> + √2(|0344>+|0455>+|0533>)
          −|1001>−|1010>+|1100>+|1154>+|1235>+|1243>−|1315>−|1332>−|1423>+|1441>−|1524>+|1552>
          +|2002>+|2020>−|2135>+|2143>+|2200>−|2254>−|2325>+|2331>−|2413>+|2442>+|2514>+|2551>
          −|3123>−|3141>+|3213>+|3242>−|3411>+|3422> + √2(|3034>+|3350>+|3505>)
          +|4115>−|4132>+|4225>+|4231>+|4311>+|4322> + √2(|4053>+|4404>+|4540>)
          −|5114>+|5151>+|5224>+|5252>+|5512>+|5521> + √2(|5045>+|5303>+|5430>)
```

In each of the three reshapings, O = 6M is, up to permutations of rows and columns, a direct sum of twelve 1 × 1
blocks (1) and twelve blocks (1/√2)[[s₁, s₂], [s₃, s₄]] with sᵢ = ±1 and s₁s₃ + s₂s₄ = 0. So O is real orthogonal
and 2-unitary, and the proof can be checked by hand. In the officers picture, each reshaping has 12 classical
officers and 24 officers in two-term Bell-type states.

## Results and how they are established

| # | Statement | Script | Status |
|---|---|---|---|
| 1 | ψ is AME(4,6): all 1- and 2-party marginals are maximally mixed | `01_state_exact.py` | exact, ℚ(√2) |
| 2 | Block structure above; 12 classical + 24 Bell-type officers in every reshaping (rows and columns); 12 product states with overlap 1/6, the maximum for any AME(4,6) state | `01` | exact |
| 3 | **Support-72 family** ψ_σ(δ, θ), δ ∈ [−1, 1], θ ∈ ℝ⁹, σ = ±1 (Eq. F72) consists of AME states. Up to diagonal local unitaries it contains every AME state with that support. It has two components, exchanged by complex conjugation, and 10 parameters modulo local unitaries | `02_families.py` | exact (Smith normal form) + numerical spot checks |
| 4 | **Support-60 face** ψ_face(θ), θ ∈ ℝ¹¹ (Eq. F60): AME, all AME states on supp ψ up to diagonal local unitaries, moduli rigid, one component; ψ = ψ_face(0) | `02` | exact + numerical |
| 5 | The families meet exactly along the 9-torus C = {θ₁₀ = θ₁₁ = 0}. **Eight branches** of the support-72 family end at ψ, four in each component, with phases ±i on the 12 vanishing entries | `02` | exact |
| 6a | At δ = −1 every branch ends at a **real** state of the opposite face (support L ∪ N) up to local phases; the simplest ends at ψ₋ = (Σ_L ε\|x⟩)/√72 + (Σ_N \|x⟩)/6 | `02` | exact |
| 6 | No state with support exactly 72 is real, even up to local phases (GF(2)). At δ = 0 the states have uniform amplitudes: √2·6M has entries in {0, ±1, ±i} with W W† = 2I, a 2-unitary complex weighing matrix. A straight path of AME states joins the uniform state with 4 imaginary entries to ψ | `02` | exact; path numerical |
| 7 | **Symmetry:** all 72 support points lie in the 360-dim −1 eigenspace of u^⊗4, u = diag(i, −1, −1, −i, −i, −i), so every state of both families (ψ included) satisfies u^⊗4\|φ⟩ = −\|φ⟩. Diagonal symmetry group of the support: Z₂ × Z₂ × Z₄ | `03_symmetry.py` | exact |
| 8 | On each face, 2³² real sign patterns give **2048 real AME states** modulo local phases (512 on the gluing locus). The support automorphism group (party permutations × level relabellings, officers to officers) has order 192, with 16 orbits on ψ's face (8 × 64, 8 × 192) and 24 on the opposite face (4 × 32, 20 × 96). No monomial map relates the two faces | `04_real_points.py` | exact |
| 9 | The real states are **isolated** modulo local orthogonal maps. Most orbits have no real first-order deformation; at the others (2 + 4 orbits), the 2-dim space of real first-order deformations is obstructed at second order | `04` | numerical (obstruction ≥ 0.35, noise ~10⁻¹⁵) |
| 10 | **Exactly 20 classes** under local unitaries and party permutations (8 + 12). Upper bound: H₁₂^⊗4 plus a monomial map merges orbits of the same face. Lower bound: the pair (N(149/8), N(3/2)), N(t) = number of eigenvalues below t of 72·JᵀJ (J = Jacobian of the six 2-party marginals), takes 20 distinct values | `05_lu_classes.py`, `07_certify_classes.py` | upper bound exact; lower bound **rigorous** (Arb interval arithmetic + Sylvester inertia) |
| 11 | First-order deformation dimensions (LU invariant): ψ and the real classes on the gluing locus: 14; the other real classes and generic face points: 11; support-72 family: 10; golden state: 4; Bruzda–Życzkowski H(0): 3. So no state here is LU-equivalent to those | `06_deformations.py` | numerical (gap 10⁻¹⁵ vs 10⁻²); for the real classes N(1/2) = 141 + dim is certified in `07` |
| 12 | Discovery: a random-start search in the 360-dim sector of u^⊗4 finds AME states within a few restarts; most hits have deformation dimension 10, like the family | `08_sector_search.py` | stochastic demonstration |

**Real classes** (data/real_classes.json). Face +: representative = ψ with the listed coefficients 1/6 negated.
Face −: representative = ψ₋ with the listed coefficients 1/6 negated.

| class | face | flipped | states | automorphism orbits | on C | deformation dim | N(149/8) | N(3/2) |
|---|---|---|---|---|---|---|---|---|
| 1 | + | – | 128 | 64 + 64 | yes | 14 | 1048 | 157 |
| 2 | + | 0344, 3505 | 384 | 192 + 192 | yes | 14 | 1046 | 157 |
| 3 | + | 5303 | 128 | 64 + 64 | no | 11 | 1050 | 153 |
| 4 | + | 5430 | 128 | 64 + 64 | no | 11 | 1044 | 152 |
| 5 | + | 5045 | 128 | 64 + 64 | no | 11 | 1053 | 153 |
| 6 | + | 0533 | 384 | 192 + 192 | no | 11 | 1055 | 153 |
| 7 | + | 0455 | 384 | 192 + 192 | no | 11 | 1052 | 153 |
| 8 | + | 0344 | 384 | 192 + 192 | no | 11 | 1052 | 152 |
| 9 | - | – | 288 | 96 + 96 + 96 | yes | 14 | 1021 | 157 |
| 10 | - | 0434, 3055 | 96 | 32 + 32 + 32 | yes | 14 | 1029 | 173 |
| 11 | - | 0434, 5340 | 96 | 96 | yes | 14 | 1025 | 159 |
| 12 | - | 0353, 4044 | 32 | 32 | yes | 14 | 985 | 155 |
| 13 | - | 0353 | 288 | 96 + 96 + 96 | no | 11 | 1004 | 152 |
| 14 | - | 0434 | 288 | 96 + 96 + 96 | no | 11 | 1032 | 160 |
| 15 | - | 0545 | 288 | 96 + 96 + 96 | no | 11 | 1022 | 152 |
| 16 | - | 3055 | 288 | 96 + 96 + 96 | no | 11 | 1026 | 152 |
| 17 | - | 0434, 5405 | 96 | 96 | no | 11 | 1038 | 168 |
| 18 | - | 5033 | 96 | 96 | no | 11 | 1020 | 152 |
| 19 | - | 5340 | 96 | 96 | no | 11 | 1012 | 152 |
| 20 | - | 5405 | 96 | 96 | no | 11 | 1030 | 152 |

## The closed forms

With the sets L (48 states, coefficient ε = ±1 in ψ), O (the 12 states with coefficient 1/6) and N (12 further states),
and the integer phase matrix K (data/family.json):

```
ψ_σ(δ, θ) = (1/√72) [ Σ_L ε e^{i(Kθ)} |x>  +  √(1+δ) Σ_O e^{i(Kθ)} |x>  +  iσ √(1−δ) Σ_N e^{i(Kθ)} |x> ]
ψ_face(θ) = (1/√72) Σ_L ε e^{i(K'θ)} |x>  +  (1/6) Σ_O e^{i(K'θ)} |x>
```

Here K′ is K restricted to L ∪ O plus two columns: θ₁₀ acts on |3505⟩ alone and θ₁₁ on |5045⟩ alone. Then
ψ = ψ_face(0) = ψ_±(1, 0). On the support of ψ every reshaping splits into 1 × 1 and 2 × 2 blocks, so the AME
conditions reduce to linear equations for the moduli and congruences φ_a − φ_b − φ_c + φ_e ≡ π (mod 2π) for the phases.
Real solutions are the solutions of the same congruences over GF(2).

`wolfram/verify_families.wls` checks symbolically, in Mathematica, that both closed forms are AME for all parameter
values (δ = cos 2α).

## How the state was found

The input was only N = 4 and d = 6. We searched for AME states inside eigenspaces of local symmetries
U^⊗4 with U = V diag(ω^e) V†, where V is Haar random and the finite-order spectra e pass simple necessary
conditions (121 classes). One sector, with spectrum {i, −1, −1, −i, −i, −i} on every party and eigenvalue −1
(dimension 360), produced an AME state. In the eigenbasis of the symmetry it has support 348. Minimising the support
over local frames inside the commutant U(1) × U(2) × U(3) of u (and over U(6)) gave support 72, an orthogonal array
of strength 2 and index 2. On it the AME conditions are linear, which gives the family in closed form. At δ = ±1
twelve entries vanish, and the remaining congruences have a solution over GF(2): ψ. Script 08 reproduces the first
step.

## Running

Python ≥ 3.10 with `numpy` and `scipy`; `07_certify_classes.py` also needs `python-flint` (Arb).

```bash
pip install -r requirements.txt
./run_all.sh            # runs 01–07 and writes results/*.log (about 6 minutes); 08 is optional
python scripts/08_sector_search.py 30 0
```

Scripts 01–07 print `[PASS]/[FAIL]` lines and exit non-zero on failure. `results/` holds the logs of a reference
run.

## Layout

```
data/psi_terms.json             exact coefficients of √72 ψ ("+1", "-1", "+sqrt2")
data/psi.npy                    ψ as a float64 vector (index = 216a + 36b + 6c + d)
data/family.json                the 72 points (type L/O/N, sign ε, row of K), face columns, branch table, θ′
data/uniform_min_imag_terms.json  a uniform-amplitude member with only 4 imaginary entries (i^k / √72)
data/real_classes.json          the 20 classes of real states (both faces)
data/reference/                 golden 2-unitary of Rather et al. (numerical, support 112);
                                matrix B of Bruzda–Życzkowski, H(0) = exp(iπB/3)/36
ame46real/                      core (cuts, AME test, blocks), lattice (Smith form, GF(2)), family (closed forms),
                                jacobian (deformations, LU-invariant spectrum), realpoints (real points, automorphisms of both faces)
scripts/01..08                  the verifications listed above
wolfram/verify_families.wls     symbolic check of the closed forms
```

## Licence

Code (`ame46real/`, `scripts/`, `wolfram/`, `run_all.sh`): MIT, see `LICENSE`. Data and results: CC BY 4.0,
see `LICENSE-DATA.md` (which also lists the sources of `data/reference/`). The manuscript in `paper/` is not covered.

## References

- D. Goyeneche et al., PRA 92, 032316 (2015): AME states and multi-unitary matrices.
- S. A. Rather et al., PRL 128, 080507 (2022): the golden AME(4,6) state; conjecture in the Supplemental Material.
- G. Rajchel-Mieldzioć et al., Rep. Prog. Phys. 89, 057601 (2026), arXiv:2508.04777: review, open problem (T2).
- W. Bruzda, K. Życzkowski, arXiv:2401.01671: two-unitary complex Hadamard matrices of order 36 (H(0)).
