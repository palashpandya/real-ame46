# A real AME(4,6) state, its families, and the classification of its real points

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
* the **classification** of all real AME states on the support of ψ: 2048 states (up to local phases) in
  **exactly 8 classes** under local unitaries and party permutations, with a rigorous (interval-arithmetic) proof
  that the 8 classes are pairwise inequivalent.

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
| 6 | No state with support exactly 72 is real, even up to local phases (GF(2)). At δ = 0 the states have uniform amplitudes: √2·6M has entries in {0, ±1, ±i} with W W† = 2I, a 2-unitary complex weighing matrix. A straight path of AME states joins the uniform state with 4 imaginary entries to ψ | `02` | exact; path numerical |
| 7 | **Symmetry:** all 72 support points lie in the 360-dim −1 eigenspace of u^⊗4, u = diag(i, −1, −1, −i, −i, −i), so every state of both families (ψ included) satisfies u^⊗4\|φ⟩ = −\|φ⟩. Diagonal symmetry group of the support: Z₂ × Z₂ × Z₄ | `03_symmetry.py` | exact |
| 8 | 2³² real sign patterns on supp ψ give **2048 real AME states** modulo local phases (512 on C). The support automorphism group (party permutations × level relabellings) has order 192 and 16 orbits (8 × 64, 8 × 192) | `04_real_points.py` | exact |
| 9 | The real states are **isolated** modulo local orthogonal maps. 14 orbits have no real first-order deformation; at the other 2, the 2-dim space of real first-order deformations is obstructed at second order | `04` | numerical (obstruction ≥ 0.35, noise ~10⁻¹⁵) |
| 10 | **Exactly 8 classes** under local unitaries and party permutations. Upper bound: H₁₂^⊗4 plus a monomial map joins the paired orbits. Lower bound: the count N(41/8) of eigenvalues below 41/8 of 72·JᵀJ (J = Jacobian of the six 2-party marginals) takes 8 distinct values: 261, 265, 255, 256, 264, 254, 257, 258 | `05_lu_classes.py`, `07_certify_classes.py` | upper bound exact; lower bound **rigorous** (Arb interval arithmetic + Sylvester inertia) |
| 11 | First-order deformation dimensions (LU invariant): ψ and class 2: 14; classes 3–8 and generic face points: 11; support-72 family: 10; golden state: 4; Bruzda–Życzkowski H(0): 3. So no state here is LU-equivalent to those | `06_deformations.py` | numerical (gap 10⁻¹⁵ vs 10⁻²); for the real classes N(1/2) = 141 + dim is certified in `07` |
| 12 | Discovery: a random-start search in the 360-dim sector of u^⊗4 finds AME states within a few restarts; most hits have deformation dimension 10, like the family | `08_sector_search.py` | stochastic demonstration |

**Real classes** (representative = ψ with the listed coefficients 1/6 negated; data/real_classes.json):

| class | flipped | states | Aut(S60) orbits | on C | deformation dim | N(41/8) |
|---|---|---|---|---|---|---|
| 1 | – (ψ) | 128 | 64 + 64 | yes | 14 | 261 |
| 2 | 0344, 3505 | 384 | 192 + 192 | yes | 14 | 265 |
| 3 | 5303 | 128 | 64 + 64 | no | 11 | 255 |
| 4 | 5430 | 128 | 64 + 64 | no | 11 | 256 |
| 5 | 5045 | 128 | 64 + 64 | no | 11 | 264 |
| 6 | 0533 | 384 | 192 + 192 | no | 11 | 254 |
| 7 | 0455 | 384 | 192 + 192 | no | 11 | 257 |
| 8 | 0344 | 384 | 192 + 192 | no | 11 | 258 |

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
./run_all.sh            # runs 01–07 and writes results/*.log (about 5 minutes); 08 is optional
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
data/real_classes.json          the 8 classes of real states
data/reference/                 golden 2-unitary of Rather et al. (numerical, support 112);
                                matrix B of Bruzda–Życzkowski, H(0) = exp(iπB/3)/36
ame46real/                      core (cuts, AME test, blocks), lattice (Smith form, GF(2)), family (closed forms),
                                jacobian (deformations, LU-invariant spectrum), realpoints (real points, Aut(S60))
scripts/01..08                  the verifications listed above
wolfram/verify_families.wls     symbolic check of the closed forms
```

## References

- D. Goyeneche et al., PRA 92, 032316 (2015): AME states and multi-unitary matrices.
- S. A. Rather et al., PRL 128, 080507 (2022): the golden AME(4,6) state; conjecture in the Supplemental Material.
- G. Rajchel-Mieldzioć et al., Rep. Prog. Phys. 89, 057601 (2026), arXiv:2508.04777: review, open problem (T2).
- W. Bruzda, K. Życzkowski, arXiv:2401.01671: two-unitary complex Hadamard matrices of order 36 (H(0)).
