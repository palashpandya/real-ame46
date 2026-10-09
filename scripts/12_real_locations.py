"""Where the 4096 real states lie in the two families (exact).

Each face is a phase torus.  Take the base real state as origin: psi on face "+" (L u O), psi_minus on face "-"
(L u N; the paper's O').  In turns theta' = theta / 2 pi, the face states up to diagonal local unitaries are
t = t0 + Kf theta' (mod 1, mod gauge), and the torus is R^11 / Lambda with
    Lambda = { theta' : Kf theta' in span_R(gauge) + Z^60 } = { theta' : B theta' in Z^m },  B = Prest Kf.
Checks:
  (1) Lambda contains Z^11; its invariant factors are printed (Lambda / Z^11).
  (2) A face state is real up to diagonal local unitaries iff 2 theta' is in Lambda, so the real states are exactly
      the 2^11 = 2048 points of order <= 2 of the torus.  This is checked as a bijection with the GF(2) count of 04.
  (3) The gluing locus C (limits of the support-72 family) is the 9-dimensional subtorus spanned by the 9 columns of
      K.  Its points of order <= 2 are exactly the 512 real states on C.  The same holds on face "-" with origin
      psi_minus.
  (4) The eight straight branches delta -> psi_sigma(delta, theta) through a real point of C end, at delta = -1, at
      real points of C' on face "-".  The bipartite graph of real points of C and C' joined by branches is printed.
"""
import itertools
from fractions import Fraction

import numpy as np

from _common import check, done
from ame46real import core, family as fam
from ame46real import realpoints as rp
from ame46real.lattice import snf


def lattice_basis(B):
    """Basis (columns, Fractions) of {y in R^k : B y in Z^m} for an integer matrix B of full column rank, and its
    invariant factors relative to Z^k."""
    D, P, Q = snf(np.asarray(B).tolist())
    k = len(B[0]); dg = [D[i][i] for i in range(k)]
    assert all(dg), "B must have full column rank"
    return [[Fraction(Q[a][i], dg[i]) for a in range(k)] for i in range(k)], sorted(abs(v) for v in dg)


def face_columns(F):
    """K restricted to the face (9 columns) plus single-officer columns completing the face torus (11 columns)."""
    G = fam.gauge_matrix(F.qs).astype(float)
    K9 = fam.K[F.qs].astype(np.int64)
    cols = [K9[:, j] for j in range(9)]
    for i in np.flatnonzero(F.officer):
        e = np.zeros(F.n, dtype=np.int64); e[i] = 1
        if np.linalg.matrix_rank(np.column_stack(cols + [e] + [G])) > np.linalg.matrix_rank(np.column_stack(cols + [G])):
            cols.append(e)
    return K9, np.column_stack(cols)


def order2_points(basis):
    k = len(basis)
    for c in itertools.product((0, 1), repeat=k):
        yield [sum(Fraction(c[i], 2) * basis[i][a] for i in range(k)) for a in range(k)]


def key_of_turns(F, t):
    y = [sum(int(p) * v for p, v in zip(row, t)) for row in F.Prest]
    return tuple(int(2 * v) % 2 for v in y)


faces = {s: rp.Face(s) for s in "+-"}
info = {}
for s, F in faces.items():
    signs, _, _ = F.real_points()
    keys = {F.key(x): i for i, x in enumerate(signs)}
    onC = {F.key(x) for x in signs if F.on_C(x)}
    K9, Kf = face_columns(F)
    t0 = [Fraction(int(v), 2) for v in F.X0]
    basis, inv = lattice_basis(F.Prest @ Kf)
    check(f"face {s}: (1) the face torus is R^11 / Lambda, Lambda contains Z^11", Kf.shape[1] == 11,
          f"Lambda / Z^11 has invariant factors {[v for v in inv if v > 1]}")
    pts = list(order2_points(basis))
    ks = []
    for th in pts:
        t = [(t0[i] + sum(int(Kf[i, a]) * th[a] for a in range(11))) % 1 for i in range(F.n)]
        assert F.real_up_to_phases(t)
        ks.append(key_of_turns(F, t))
    check(f"face {s}: (2) the 2^11 points of order <= 2 are exactly the {len(signs)} real states",
          len(set(ks)) == 2048 and set(ks) == set(keys))
    basis9, inv9 = lattice_basis(F.Prest @ K9)
    pts9 = list(order2_points(basis9)); ks9 = []
    for th in pts9:
        t = [(t0[i] + sum(int(K9[i, a]) * th[a] for a in range(9))) % 1 for i in range(F.n)]
        ks9.append(key_of_turns(F, t))
    check(f"face {s}: (3) the 2^9 points of order <= 2 of the gluing locus C are exactly the 512 real states on C",
          len(set(ks9)) == 512 and set(ks9) == onC)
    info[s] = dict(keys=keys, onC=onC, basis9=basis9, t0=t0, K9=K9)

# (4) branches from the real points of C (face +) to face -
Fp, Fm = faces["+"], faces["-"]
qm = Fm.qs
ends = {}
for th in order2_points(info["+"]["basis9"]):
    start = key_of_turns(Fp, [(info["+"]["t0"][i] + sum(int(info["+"]["K9"][i, a]) * th[a] for a in range(9))) % 1
                              for i in range(Fp.n)])
    for b in fam.F["branches"]:
        top = [2 * th[a] + Fraction(b["theta_over_pi"][a]) for a in range(9)]       # theta / pi
        t72 = fam.turns72(b["sigma"], top)
        tp = [t72[q] for q in Fp.qs]
        assert key_of_turns(Fp, tp) == start and Fp.real_up_to_phases(tp)         # delta = +1 end is the start
        tm = [t72[q] for q in qm]
        if not Fm.real_up_to_phases(tm):
            ends.setdefault(start, []).append(None); continue
        ends.setdefault(start, []).append(key_of_turns(Fm, tm))
# class labels (as in 05_lu_classes.py: automorphism orbits merged by H12^(x)4, named by the class table)
H = np.eye(6); H[1:3, 1:3] = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
table = core.load_json("real_classes.json")["classes"]
cls_of = {}
for s, F in faces.items():
    signs, _, _ = F.real_points()
    orbs = F.orbits(signs, F.automorphisms())
    orbit_of = {F.key(signs[i]): k for k, o in enumerate(orbs) for i in o}
    parent = list(range(len(orbs)))
    def find(a):
        while parent[a] != a: a = parent[a]
        return a
    for k, o in enumerate(orbs):
        w = np.einsum("ai,bj,ck,dl,ijkl->abcd", H, H, H, H, F.state(signs[o[0]]).reshape((6,) * 4)).ravel()
        x = F.to_face_signs(w)
        if x is not None:
            a, b = find(k), find(orbit_of[F.key(x)])
            if a != b: parent[a] = b
    root_cls = {find(orbit_of[F.key(F.flips(r["flip"]))]): r["cls"] for r in table if r["face"] == s}
    for kk, k in orbit_of.items():
        cls_of[(s, kk)] = root_cls[find(k)]
import collections
real_starts = collections.Counter(cls_of[("+", st)] for st, v in ends.items() if all(e is not None for e in v))
dead_starts = collections.Counter(cls_of[("+", st)] for st, v in ends.items() if all(e is None for e in v))
mixed = sum(1 for v in ends.values() if any(e is None for e in v) and any(e is not None for e in v))
real_ends = collections.Counter(cls_of[("-", e)] for e in {e for v in ends.values() for e in v if e is not None})
print(f"    starts on C whose 8 branches all end at real states: {dict(real_starts)} (class: count)")
print(f"    starts on C whose branches all end at non-real states: {dict(dead_starts)}; mixed starts: {mixed}")
print(f"    real branch ends on C': {dict(real_ends)} (class: count)")
check("(4) branch ends are either all real or all non-real; the real ones lie on C'",
      mixed == 0 and all(e in info["-"]["onC"] for v in ends.values() for e in v if e is not None))
all_real = all(e is not None for v in ends.values() for e in v)
all_onC = all(e in info["-"]["onC"] for v in ends.values() for e in v if e is not None)
deg_start = sorted({len(set(v)) for v in ends.values()})
inc = {}
for st, v in ends.items():
    for e in set(v):
        inc.setdefault(e, set()).add(st)
deg_end = sorted({len(v) for v in inc.values()})
print(f"    branch graph: {len(ends)} real points of C, {len(inc)} real points of C' reached; distinct ends per start "
      f"{deg_start}; starts per end {deg_end}")
done()
