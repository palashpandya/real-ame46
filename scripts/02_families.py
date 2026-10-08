"""The two continuous families through psi (closed forms, Eqs. F72 and F60) and how they meet.

Exact (integer / rational arithmetic):
  (a) support 72 = L u O u N: every cut is a direct sum of 18 2x2 blocks; the moduli equations have a 1-dimensional
      affine solution set (delta); the 54 phase congruences have a 30-dim real kernel = span(K, local phases),
      i.e. 9 phase parameters mod LU; Smith form has one invariant factor 2 -> exactly 2 components, exchanged by
      complex conjugation, labelled by sigma = +-1.
  (b) support 60 = L u O (the face): moduli rigid; 36 congruences, kernel 32 = span(K', local phases) -> 11
      parameters; all invariant factors 1 -> one component.  psi = psi_face(0).
  (c) gluing: a face state is a limit of support-72 AME states iff theta_10, theta_11 = 0 (mod 2 pi).
  (d) eight branches: with psi's 60 phases fixed, the 54 x 12 system for the vanishing phases has invariant factors
      1^9 2^3 -> exactly 8 solutions, all +-i; they are the 8 rows of the branch table (4 per component).
  (e) no real state has support exactly L u O u N, even up to local phases (GF(2)); psi's face does have real points.
  (g) the other end: at delta = -1 every branch ends at a REAL state of the opposite face (support L u N) up to
      local phases; the simplest one ends at psi_minus = (sum_L eps|x>)/sqrt72 + (sum_N |x>)/6;
  (f) uniform amplitudes: at delta = 0, sqrt2 * 6M has entries in {0, +-1, +-i} and W W^dag = 2 I in all three cuts
      (a 2-unitary complex weighing matrix); the same for the stored state with only 4 imaginary entries, which
      equals psi_sigma(0, theta') up to local phases.
Numerical: AME of both families at random parameters, along all 8 branches, and along the path theta' -> psi.
"""
from fractions import Fraction as Fr

import numpy as np

from _common import check, done
from ame46real import core, family as fam
from ame46real.realpoints import Face
from ame46real.lattice import all_solutions_mod1, invariant_factors, snf, solve_mod1, gf2_solvable

psi = core.psi()
t60 = [Fr(0) if fam.EPS[q] > 0 else Fr(1, 2) for q in fam.S60]          # phases of psi (turns) on S60

# (a) support 72 ---------------------------------------------------------------------------------------------
M72, quads72, singles72 = fam.congruence_matrix(fam.S72)
G72 = fam.gauge_matrix(fam.S72)
check("(a) S72: each cut is 18 blocks 2x2", len(quads72) == 54 and not singles72)
rows = []
for a, b, c, e in quads72:
    for u, w, sgn, rhs in ((a, e, -1, 0), (b, c, -1, 0), (a, b, 1, 2)):
        row = np.zeros(72); row[fam.POS[u]] = 1; row[fam.POS[w]] = sgn; rows.append(np.append(row, rhs))
A = np.array(rows)
s_delta = np.array([0 if t == "L" else (1 if t == "O" else -1) for t in fam.TYPE])
check("(a) moduli s = 72|a|^2: affine solution set is {1 + delta (O), 1 - delta (N), 1 (L)}",
      72 - np.linalg.matrix_rank(A[:, :72]) == 1 and np.allclose(A[:, :72] @ np.ones(72), A[:, 72])
      and np.allclose(A[:, :72] @ s_delta, 0))
rk = np.linalg.matrix_rank(M72); rG = np.linalg.matrix_rank(G72)
check("(a) phase congruences: kernel 30, local phases 21 -> 9 parameters; K spans the rest",
      72 - rk == 30 and rG == 21 and not np.any(M72 @ fam.K) and np.linalg.matrix_rank(np.hstack([fam.K, G72])) == 30,
      f"rank M = {rk}")
comps = fam.Components(M72)
check("(a) invariant factors -> exactly 2 components", comps.count == 2, f"factors {sorted(set(comps.d))}")
bp, bm = fam.base_turns(1), fam.base_turns(-1)
for s, b in ((1, bp), (-1, bm)):
    Mb = [sum(int(M72[r, q]) * b[q] for q in range(72)) % 1 for r in range(len(quads72))]
    check(f"(a) sigma = {s:+d} base point solves the congruences", all(v == Fr(1, 2) for v in Mb))
check("(a) sigma = +1 and sigma = -1 lie in different components", not comps.same(bp, bm))
check("(a) complex conjugation maps component sigma=+1 to sigma=-1", comps.same([-v % 1 for v in bp], bm))

# (b) the face ------------------------------------------------------------------------------------------------
M60, quads60, singles60 = fam.congruence_matrix(fam.S60)
G60 = fam.gauge_matrix(fam.S60); Kf = fam.face_K()
check("(b) S60: per cut 12 blocks 2x2 (L) + 12 blocks 1x1 (O)", len(quads60) == 36 and len(singles60) == 36)
check("(b) moduli rigid: officers alone in row and column (s = 2), L blocks force s = 1",
      all(fam.TYPE[fam.POS[x]] == "O" for x in singles60) and all(fam.TYPE[fam.POS[x]] == "L" for q in quads60 for x in q))
rk60 = np.linalg.matrix_rank(M60)
check("(b) kernel 32, local phases 21 -> 11 parameters; K' spans; one component",
      60 - rk60 == 32 and np.linalg.matrix_rank(G60) == 21 and not np.any(M60 @ Kf)
      and np.linalg.matrix_rank(np.hstack([Kf, G60])) == 32 and set(invariant_factors(M60.tolist())) == {1})
check("(b) psi = psi_face(0) = psi_{+}(1, 0) = psi_{-}(1, 0)",
      np.allclose(fam.psi_face(np.zeros(11)), psi) and np.allclose(fam.psi72(1, 1.0, np.zeros(9)), psi)
      and np.allclose(fam.psi72(-1, 1.0, np.zeros(9)), psi))

# (c) gluing locus --------------------------------------------------------------------------------------------
Au = M72[:, fam.NPTS]                                       # unknown phases of the vanishing points
D_, P_, _ = snf(Au.tolist())
P_ = np.array(P_, dtype=object)
R = (P_ @ (M72[:, fam.S60].astype(object) @ Kf.astype(object)))[12:]   # theta-dependence of the solvability conditions
check("(c) solvability for the N phases depends only on theta_10, theta_11", not np.any(R[:, :9]))
check("(c) ... and holds iff theta_10, theta_11 in 2 pi Z (invariant factors 1, 1)",
      invariant_factors(R[:, 9:].astype(np.int64).tolist()) == [1, 1])

# (d) the eight branches through psi --------------------------------------------------------------------------
rhs = [Fr(1, 2) - sum(int(M72[r, q]) * t for q, t in zip(fam.S60, t60)) for r in range(len(quads72))]
sols = all_solutions_mod1(Au.tolist(), rhs)
check("(d) 54 x 12 system: invariant factors 1^9 2^3 -> 8 solutions, all phases +-i",
      sorted(invariant_factors(Au.tolist())) == [1] * 9 + [2] * 3 and len(sols) == 8
      and all(v in (Fr(1, 4), Fr(3, 4)) for y in sols for v in y))
officer = np.array([fam.TYPE[q] == "O" for q in fam.S60])
errs = []
for y in sols:
    for delta in (0.999, 0.5, 0.0, -0.5, -1.0):
        v = np.zeros(core.DIM, complex)
        v[[fam.FLAT[q] for q in fam.S60]] = psi[[fam.FLAT[q] for q in fam.S60]] * np.where(officer, np.sqrt((1 + delta) / 2), 1)
        v[[fam.FLAT[q] for q in fam.NPTS]] = np.exp(2j * np.pi * np.array([float(t) for t in y])) * np.sqrt((1 - delta) / 72)
        errs.append(core.ame_error(v))
check("(d) each solution gives an AME curve delta in [-1, 1] ending at psi (numerical)", max(errs) < 1e-14, f"max error {max(errs):.1e}")
# the table rows (sigma, theta) of data/family.json are these 8 branches
# local phases g fixing every phase of psi (G60 g in Z^60) must act trivially on N for the branches to be well defined
D6, P6, Q6 = snf(G60.tolist()); r6 = sum(1 for i in range(min(len(D6), len(D6[0]))) if D6[i][i]); Q6 = np.array(Q6)
acts = [[Fr(int(v), D6[i][i]) % 1 for v in G72[fam.NPTS] @ Q6[:, i]] for i in range(r6)]   # generators (turns on N)
check("(d) the local phases fixing psi (a finite group, factors " + str([D6[i][i] for i in range(r6) if D6[i][i] > 1])
      + ") act trivially on N", not np.any(G72[fam.NPTS] @ Q6[:, r6:]) and all(v == 0 for a in acts for v in a))
hit, per_sigma = set(), {1: 0, -1: 0}
for b in fam.F["branches"]:
    T = fam.turns72(b["sigma"], b["theta_over_pi"])
    g = solve_mod1(G60.tolist(), [(T[q] - t) % 1 for q, t in zip(fam.S60, t60)])
    if g is None:
        continue
    y = tuple((T[q] - sum(int(G72[q, k]) * g[k] for k in range(24))) % 1 for q in fam.NPTS)
    if list(y) in sols:
        hit.add(y); per_sigma[b["sigma"]] += 1
check("(d) the 8 table rows are the 8 branches, 4 in each component", len(hit) == 8 and per_sigma == {1: 4, -1: 4})

# (g) delta = -1 ends -------------------------------------------------------------------------------------------
Fm = Face("-")
ends = [Fm.real_up_to_phases([fam.turns72(b["sigma"], b["theta_over_pi"])[q] for q in Fm.qs]) for b in fam.F["branches"]]
check("(g) all 8 branches end (delta = -1) at real states of the opposite face, up to local phases", all(ends))
Mm = fam.congruence_matrix(Fm.qs)[0]
check("(g) the opposite face is also an 11-parameter phase torus with one component (36 congruences, kernel 32)",
      Mm.shape[0] == 36 and 60 - np.linalg.matrix_rank(Mm) == 32 and np.linalg.matrix_rank(fam.gauge_matrix(Fm.qs)) == 21
      and set(invariant_factors(Mm.tolist())) == {1})
tm = [Fr(0) if x == 0 else Fr(1, 2) for x in Fm.X0]
check("(g) psi_+(-1, 0) and psi_-(-1, 0) equal psi_minus up to local phases; psi_minus is AME",
      all(fam.equal_up_to_local_phases([fam.turns72(sg, [0] * 9)[q] for q in Fm.qs], tm, Fm.qs) for sg in (1, -1))
      and core.ame_error(Fm.state(Fm.X0)) < 1e-14)

# (e) reality ------------------------------------------------------------------------------------------------
check("(e) no real point on support 72: M x = 1 has no solution over GF(2)", not gf2_solvable(M72 % 2, np.ones(54, int)))
check("(e) real points exist on the face: psi's signs solve M60 x = 1 over GF(2)", gf2_solvable(M60 % 2, np.ones(36, int)))

# (f) uniform amplitudes -------------------------------------------------------------------------------------
def weighing(v72):
    """v72: sqrt72 * state with entries in {0, +-1, +-i}; check W W^dag = 2 I in all cuts (integer arithmetic)."""
    W = np.round(v72).astype(complex)
    if not np.allclose(W, v72) or not set(np.abs(W[np.abs(W) > 0]).round(12)) == {1.0}:
        return False
    return all(np.array_equal(core.cut_matrix(W, c) @ core.cut_matrix(W, c).conj().T, 2 * np.eye(36)) for c in core.CUTS)


for s in (1, -1):
    check(f"(f) psi_{s:+d}(0, 0): sqrt2*6M in {{0, +-1, +-i}}, W W^dag = 2 I in all cuts",
          weighing(np.sqrt(72) * fam.psi72(s, 0.0, np.zeros(9))))
uni = core.load_json("uniform_min_imag_terms.json")["terms"]
U = np.zeros(core.DIM, complex)
for x, c in uni:
    U[core.index(core.tup(x))] = {"+1": 1, "-1": -1, "+i": 1j, "-i": -1j}[c]
check("(f) stored uniform state (68 x +-1, 4 x +-i): exact weighing 2-unitary", weighing(U) and sum("i" in c for _, c in uni) == 4)
cf = fam.F["uniform_closed_form"]
tU = [Fr(int(round(np.angle(U[fam.FLAT[q]]) / (np.pi / 2))) % 4, 4) for q in fam.S72]
check("(f) it equals psi_sigma(0, theta') up to local phases (exact)",
      fam.equal_up_to_local_phases(tU, fam.turns72(cf["sigma"], cf["theta_over_pi"]), fam.S72))
th = np.pi * np.array([float(Fr(v)) for v in cf["theta_over_pi"]])
path = [core.ame_error(fam.psi72(cf["sigma"], t, (1 - t) * th)) for t in np.linspace(0, 1, 21)]
check("(f) straight path t -> psi_sigma(t, (1-t) theta') is AME and ends at psi (numerical)",
      max(path) < 1e-14 and np.allclose(fam.psi72(cf["sigma"], 1.0, np.zeros(9)), psi), f"max error {max(path):.1e}")

# numerical AME at random parameters ---------------------------------------------------------------------------
rng = np.random.default_rng(2026)
e72 = max(core.ame_error(fam.psi72(s, rng.uniform(-1, 1), rng.uniform(0, 2 * np.pi, 9))) for s in (1, -1) for _ in range(20))
e60 = max(core.ame_error(fam.psi_face(rng.uniform(0, 2 * np.pi, 11))) for _ in range(20))
check("random members of both families are AME (numerical, 40 + 20 points)", max(e72, e60) < 1e-14,
      f"max errors {e72:.1e}, {e60:.1e}")
th9 = rng.uniform(0, 2 * np.pi, 9)
check("psi_sigma(1, theta) = psi_face(theta, 0, 0)", np.allclose(fam.psi72(1, 1.0, th9), fam.psi_face(np.append(th9, [0, 0]))))
done()
