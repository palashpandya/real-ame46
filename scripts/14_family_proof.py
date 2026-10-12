"""The computational facts used in the proof of Proposition 2 (alpha/beta form of the families, paper Sec. IV).

E = 48 states of psi with coefficient +-1/sqrt72, P = 12 with 1/6, P' = the 12 further states of the 72-point support;
P and P' are split into the quartets P_1, P_2, P_3 and P'_1, P'_2, P'_3 of Eq. (quartets) (paper/inputs/family_sets.tex).

Exact (integer / rational arithmetic):
  (1) Lemma 1 applies: on E u P u P' every cut is a direct sum of 12 full 2x2 E-blocks and 6 full mixed blocks
      (P on one diagonal, P' on the other); on E u P and E u P' of 12 E-blocks and 12 1x1 blocks.
  (2) the 108 moduli equations of the 36 E-blocks have rank 48 (only solution s = 1 on E).
  (3) the 36 E-congruences in the 48 phases on E: rank 28, all invariant factors 1; the diagonal local phases on E
      have rank 20 and solve the homogeneous system.
  (4) mixed blocks: every pair within a quartet lies in exactly one of them; for every (k, l) exactly two contain
      states of P_k and P'_l, and their pairs partition P_k and P'_l; partner pairs have the same first indices.
  (5) the 18 mixed congruences in the 24 phases of P u P': rank 14, invariant factors 1^13 2.
  (6) omega: the 9 mixed blocks without first index 0 contain each x in P u P' with x_1 != 0 exactly twice and no
      other state; omega is unchanged by diagonal local phases; omega(0, +-pi/2) = +-pi/2.
  (7) the diagonal local unitaries that fix the coefficients on E act on P u P' by a common phase and by the sign
      change of P_3 u P' (and products).
  (8) quartet sums: for random alpha with equal quartet sums, Eq. (mixed) has exactly 8 solutions beta, two on each
      quartet P'_l (differing by pi on P'_l), 4 with omega = +pi/2 and 4 with -pi/2; for unequal sums none.  The same
      with alpha and beta exchanged.
  (9) no solution over GF(2) of the mixed congruences (no real state with support E u P u P').
Numerical: Psi_delta(alpha, beta) is AME at random points of both components, and Psi_1(0) = psi.
"""
import itertools
import os
import random
import re
from fractions import Fraction as Fr

import numpy as np

from _common import check, done
from ame46real import core, family as fam
from ame46real.lattice import all_solutions_mod1, gf2_solvable, invariant_factors, snf

rng = random.Random(7)
st = open(os.path.join(core.DATA, os.pardir, "paper", "inputs", "family_sets.tex")).read()
QUART = {("P" + pr + k): [core.tup(x) for x in re.findall(r"\d{4}", body)]
         for pr, k, body in re.findall(r"P('?)_(\d)&=\\\{(.*?)\\\}", st)}
PK = [QUART[f"P{k}"] for k in (1, 2, 3)]
QK = [QUART[f"P'{k}"] for k in (1, 2, 3)]
PP, QQ = sum(PK, []), sum(QK, [])
TYPE = {fam.POINTS[q]: fam.TYPE[q] for q in fam.S72}
E = [x for x in fam.POINTS if TYPE[x] == "L"]
check("quartets of the paper partition P (data: O) and P' (data: N)",
      sorted(PP) == sorted(x for x in TYPE if TYPE[x] == "O") and sorted(QQ) == sorted(x for x in TYPE if TYPE[x] == "N"))
KOF = {x: k for k in range(3) for x in PK[k]} | {x: k for k in range(3) for x in QK[k]}


# ---- (1) blocks -------------------------------------------------------------------------------------------------
def block_list(points):
    """[(cut, R, C, states)] for every connected component of the row/column graph, with its states."""
    out = []
    for cut in core.CUTS:
        for R, C in core.blocks(points, cut):
            out.append((cut, R, C, [core.point(r, c, cut) for r in R for c in C]))
    return out


S72 = E + PP + QQ
B72 = block_list(S72)
full = all(len(R) == len(C) and set(xs) <= set(S72) for _, R, C, xs in B72)
Eb = [b for b in B72 if all(TYPE[x] == "L" for x in b[3])]
Mb = [b for b in B72 if b not in Eb]
mixed_ok = all(len(b[1]) == 2 and {TYPE[b[3][0]], TYPE[b[3][3]]} in ({"O"}, {"N"})
               and {TYPE[b[3][1]], TYPE[b[3][2]]} == ({"O", "N"} - {TYPE[b[3][0]]}) for b in Mb)
check("E u P u P': each cut = 12 full 2x2 E-blocks + 6 full mixed blocks (P on one diagonal, P' on the other)",
      full and all(len(b[1]) == 2 for b in B72) and len(Eb) == 36 and len(Mb) == 18 and mixed_ok
      and all(sum(1 for b in Eb if b[0] == c) == 12 for c in core.CUTS))
for name, S in (("E u P", E + PP), ("E u P'", E + QQ)):
    B = block_list(S)
    ok = all(len(R) == len(C) and set(xs) <= set(S) for _, R, C, xs in B)
    n2 = [b for b in B if len(b[1]) == 2]; n1 = [b for b in B if len(b[1]) == 1]
    check(f"{name}: each cut = 12 full 2x2 E-blocks + 12 blocks 1x1", ok and len(n2) == 36 and len(n1) == 36
          and all(all(TYPE[x] == "L" for x in b[3]) for b in n2) and len(B) == 72)

# ---- (2) moduli on E --------------------------------------------------------------------------------------------
iE = {x: i for i, x in enumerate(E)}
rows = []
for _, _, _, (a, b, c, e) in Eb:                                 # [[a, b], [c, e]]
    for u, v, s in ((a, e, -1), (b, c, -1), (a, b, 1)):
        r = np.zeros(48, dtype=np.int64); r[iE[u]] += 1; r[iE[v]] += s; rows.append(r)
check("moduli equations of the 36 E-blocks: 108 equations of rank 48 (only solution s = 1 on E)",
      len(rows) == 108 and np.linalg.matrix_rank(np.array(rows, float)) == 48)

# ---- (3) E-congruences and local phases -------------------------------------------------------------------------
ME = np.zeros((36, 48), dtype=np.int64)
for i, (_, _, _, (a, b, c, e)) in enumerate(Eb):
    ME[i, iE[a]] += 1; ME[i, iE[b]] -= 1; ME[i, iE[c]] -= 1; ME[i, iE[e]] += 1
Lmat = lambda pts: np.array([[int(x[j] == v) for j in range(core.N) for v in range(core.D)] for x in pts], dtype=np.int64)
LE = Lmat(E)
check("E-congruences: rank 28, all invariant factors 1", invariant_factors(ME.tolist()) == [1] * 28)
check("diagonal local phases on E: rank 20, and they solve the homogeneous E-congruences",
      np.linalg.matrix_rank(LE.astype(float)) == 20 and not (ME @ LE).any())

# ---- (4) structure of the mixed blocks --------------------------------------------------------------------------
pairs = []                                                       # (pair in P, pair in P')
for _, _, _, (a, b, c, e) in Mb:
    d1, d2 = frozenset((a, e)), frozenset((b, c))
    pairs.append((d1, d2) if TYPE[a] == "O" else (d2, d1))
allP = [p for p, _ in pairs]; allQ = [q for _, q in pairs]
K4 = all(sorted(map(sorted, (p for p in allP if KOF[next(iter(p))] == k))) ==
         sorted(sorted(s) for s in map(frozenset, itertools.combinations(PK[k], 2))) for k in range(3)) \
     and all(sorted(map(sorted, (q for q in allQ if KOF[next(iter(q))] == k))) ==
             sorted(sorted(s) for s in map(frozenset, itertools.combinations(QK[k], 2))) for k in range(3))
same_quartet = all(len({KOF[x] for x in p}) == 1 and len({KOF[x] for x in q}) == 1 for p, q in pairs)
check("every pair within a quartet lies in exactly one mixed block", same_quartet and K4)
kl_ok = True
for k in range(3):
    for l in range(3):
        bl = [(p, q) for p, q in pairs if KOF[next(iter(p))] == k and KOF[next(iter(q))] == l]
        kl_ok &= len(bl) == 2 and set().union(*(p for p, _ in bl)) == set(PK[k]) and set().union(*(q for _, q in bl)) == set(QK[l])
check("for every (k, l) exactly two mixed blocks contain states of P_k and P'_l; their pairs partition P_k and P'_l",
      kl_ok)
check("partner pairs have the same first indices; each quartet has one state with each first index 0, 3, 4, 5",
      all({x[0] for x in p} == {x[0] for x in q} for p, q in pairs)
      and all(sorted(x[0] for x in Qk) == [0, 3, 4, 5] for Qk in PK + QK))

# ---- (5) mixed congruences --------------------------------------------------------------------------------------
iV = {x: i for i, x in enumerate(PP + QQ)}
MM = np.zeros((18, 24), dtype=np.int64)                          # alpha_a + alpha_b - beta_a' - beta_b' = pi
for i, (p, q) in enumerate(pairs):
    for x in p: MM[i, iV[x]] += 1
    for x in q: MM[i, iV[x]] -= 1
check("mixed congruences: rank 14, invariant factors 1^13 2", invariant_factors(MM.tolist()) == [1] * 13 + [2])

# ---- (6) omega ----------------------------------------------------------------------------------------------------
nine = [i for i, (p, q) in enumerate(pairs) if all(x[0] != 0 for x in p | q)]
cover = MM[nine].sum(axis=0)
w = np.array([-1 if x in PP else 1 for x in PP + QQ]) * np.array([int(x[0] != 0) for x in PP + QQ])   # omega
check("the 9 mixed blocks without first index 0 contain each x with x_1 != 0 exactly twice, no other state",
      len(nine) == 9 and np.array_equal(-cover, 2 * w))
check("omega is unchanged by diagonal local phases", not (w @ Lmat(PP + QQ)).any())
om = lambda gam: sum(int(c) * Fr(g) for c, g in zip(w, gam)) % 1                  # phases in turns
check("omega(Psi(0, +pi/2)) = +pi/2 and omega(Psi(0, -pi/2)) = -pi/2",
      om([0] * 12 + [Fr(1, 4)] * 12) == Fr(1, 4) and om([0] * 12 + [Fr(3, 4)] * 12) == Fr(3, 4))

# ---- (7) residual diagonal gauge ----------------------------------------------------------------------------------
Dm, _, Q = snf(LE.tolist())
Q = np.array(Q, dtype=object)
d = [Dm[i][i] for i in range(min(len(Dm), len(Dm[0]))) if Dm[i][i]]
LPQ = Lmat(PP + QQ)
cont = [LPQ @ Q[:, i] for i in range(len(d), 24)]                # continuous part of ker(L_E: T^24 -> T^48)
fin = [[Fr(int(v), d[i]) % 1 for v in LPQ @ Q[:, i]] for i in range(len(d)) if d[i] > 1]
flip = [Fr(0)] * 8 + [Fr(1, 2)] * 16                             # pi on P_3 u P'
circ = lambda v: len(set(v)) == 1
check("E-fixing diagonal gauge on P u P': circle of common phases and the sign change of P_3 u P'",
      d == [1] * 16 + [2] * 4 and np.linalg.matrix_rank(np.array(cont, float)) == 1 and all(circ(v) for v in cont)
      and all(circ(v) or circ([(a - b) % 1 for a, b in zip(v, flip)]) for v in fin)
      and any(circ([(a - b) % 1 for a, b in zip(v, flip)]) for v in fin))


# ---- (8) quartet sums ------------------------------------------------------------------------------------------
Ma, Mq = MM[:, :12], -MM[:, 12:]                                 # Ma alpha - Mq beta = 1/2
def betas(alpha):
    return all_solutions_mod1(Mq.tolist(), [(sum(int(Ma[i, j]) * alpha[j] for j in range(12)) - Fr(1, 2)) % 1 for i in range(18)])
def alphas(beta):
    return all_solutions_mod1(Ma.tolist(), [(sum(int(Mq[i, j]) * beta[j] for j in range(12)) + Fr(1, 2)) % 1 for i in range(18)])


def rand_side(equal):
    v = [Fr(rng.randrange(48), 48) for _ in range(12)]
    s = [sum(v[4 * k:4 * k + 4]) % 1 for k in range(3)]
    if equal:
        v[7] = (v[7] + s[0] - s[1]) % 1; v[11] = (v[11] + s[0] - s[2]) % 1
    elif len(set(s)) == 1:
        v[7] = (v[7] + Fr(1, 48)) % 1
    return v


ok8 = True
for _ in range(20):
    a = rand_side(True); bs = betas(a)
    pairs_ok = all(any(all((b2[j] - b1[j]) % 1 == (Fr(1, 2) if j // 4 == l else 0) for j in range(12)) for b2 in bs)
                   for b1 in bs for l in range(3))
    oms = sorted(om(list(a) + list(b)) for b in bs)
    ok8 &= len(bs) == 8 and pairs_ok and oms == [Fr(1, 4)] * 4 + [Fr(3, 4)] * 4
    ok8 &= betas(rand_side(False)) == []
check("equal quartet sums of alpha: exactly 8 solutions beta (2 per quartet, differing by pi), 4 in each component; "
      "unequal sums: none", ok8)
ok8 = True
for _ in range(20):
    b = rand_side(True); As = alphas(b)
    ok8 &= len(As) == 8 and sorted(om(list(a) + list(b)) for a in As) == [Fr(1, 4)] * 4 + [Fr(3, 4)] * 4
    ok8 &= alphas(rand_side(False)) == []
check("the same with alpha and beta exchanged", ok8)

# ---- (9) GF(2) ------------------------------------------------------------------------------------------------------
check("the mixed congruences have no solution over GF(2)", not gf2_solvable(MM % 2, np.ones(18, dtype=np.int64)))


# ---- numerical -------------------------------------------------------------------------------------------------------
def Psi(delta, alpha, beta):
    v = np.zeros(core.DIM, dtype=complex)
    eps = {fam.POINTS[q]: fam.EPS[q] for q in fam.S72}
    for x in E: v[core.index(x)] = eps[x]
    for x, t in zip(PP, alpha): v[core.index(x)] = np.sqrt(1 + delta) * np.exp(2j * np.pi * float(t))
    for x, t in zip(QQ, beta): v[core.index(x)] = np.sqrt(1 - delta) * np.exp(2j * np.pi * float(t))
    return v / np.sqrt(72)


check("Psi_1(0) = psi", np.allclose(Psi(1, [0] * 12, [0] * 12), core.psi()))
errs = []
for _ in range(10):
    a = rand_side(True)
    for b in betas(a):
        errs.append(core.ame_error(Psi(rng.uniform(-0.99, 0.99), a, b)))
check("Psi_delta(alpha, beta) is AME at random points of both components", max(errs) < 1e-12, f"max error {max(errs):.1e}")
done()
