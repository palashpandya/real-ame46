"""Real points of the support-60 face and the action of the support automorphism group Aut(S60).

A real state with support exactly S60 = L u O has signs x in GF(2)^60 (1 = minus) with Mg x = 1, where Mg is the
block-congruence matrix mod 2.  Two sign patterns give the same state up to diagonal local phases iff
Pc (x - x') = 0 (mod 2), where Pc are the rows of the Smith transform P of the gauge matrix beyond its rank.
"""
import collections
import itertools

import numpy as np

from .core import D, N, DIM
from . import family as fam
from .lattice import snf, gf2_nullspace, gf2_rank

X = [fam.POINTS[q] for q in fam.S60]
POS = {x: i for i, x in enumerate(X)}
FLAT = [fam.FLAT[q] for q in fam.S60]
MOD = np.array([1.0 if fam.TYPE[q] == "L" else 2.0 for q in fam.S60])     # 72 |amplitude|^2
X0 = (fam.EPS[fam.S60] < 0).astype(np.int64)                             # signs of psi


def state(x):
    """Real face state with sign pattern x (on S60)."""
    v = np.zeros(DIM)
    v[FLAT] = np.sqrt(MOD / 72) * (1 - 2 * np.asarray(x))
    return v


def _cokernel_rows(B):
    D_, P, _ = snf(np.asarray(B).tolist())
    r = sum(1 for i in range(min(len(D_), len(D_[0]))) if D_[i][i])
    return np.array(P, dtype=np.int64)[r:] % 2


Mg = fam.congruence_matrix(fam.S60)[0] % 2
G60 = fam.gauge_matrix(fam.S60)
Pc = _cokernel_rows(G60)                                                  # key: real point modulo local phases
P2c = _cokernel_rows(np.hstack([fam.face_K()[:, :9], G60]))              # (x - x0) on the gluing locus C


def key(x):
    return tuple((Pc @ np.asarray(x)) % 2)


def real_points():
    """Returns (sign patterns of the distinct real points, dim of the GF(2) solution space, dim of the gauge part)."""
    assert not np.any((Mg @ X0 - 1) % 2)
    ker = gf2_nullspace(Mg)
    comp, img = [], []
    for v in ker:                                    # complement of the gauge part inside ker(Mg)
        trial = img + [(Pc @ v) % 2]
        if gf2_rank(np.array(trial)) == len(trial):
            comp.append(v); img = trial
    signs = [(X0 + sum(c * v for c, v in zip(cs, comp))) % 2 for cs in itertools.product((0, 1), repeat=len(comp))]
    return signs, len(ker), len(ker) - len(comp)


def on_C(x):
    return not np.any((P2c @ ((np.asarray(x) - X0) % 2)) % 2)


def theta_points():
    """Distinct real points reached by theta_j in {0, pi} in the face closed form."""
    Kf = fam.face_K()
    return {key((X0 + Kf @ np.array(t)) % 2) for t in itertools.product((0, 1), repeat=11)}


# ---------- support automorphisms: party permutations x level relabellings mapping S60 onto itself ----------
def relabelings(Y, target, sp):
    """All level relabellings pi (one per party) with {pi(Y[:, sp])} = target, found by backtracking."""
    Y = np.asarray(Y)[:, list(sp)]
    proj = {sub: {tuple(t[j] for j in sub) for t in target} for r in range(2, N + 1)
            for sub in itertools.combinations(range(N), r)}
    out, pi, used = [], [dict() for _ in range(N)], [set() for _ in range(N)]
    order = [(j, v) for v in range(D) for j in range(N)]

    def ok():
        for y in Y:
            s = tuple(j for j in range(N) if y[j] in pi[j])
            if len(s) >= 2 and tuple(pi[j][y[j]] for j in s) not in proj[s]:
                return False
        return True

    def rec(i):
        if i == len(order):
            if {tuple(pi[j][y[j]] for j in range(N)) for y in Y} == target:
                out.append([dict(p) for p in pi])
            return
        j, v = order[i]
        for w in range(D):
            if w not in used[j]:
                pi[j][v] = w; used[j].add(w)
                if ok():
                    rec(i + 1)
                del pi[j][v]; used[j].discard(w)
    rec(0)
    return out


def automorphisms():
    """Aut(S60) as permutations of the 60 points (img[i] = image of point i)."""
    S = set(X)
    perms = []
    for sp in itertools.permutations(range(N)):
        for pi in relabelings(X, S, sp):
            img = np.array([POS[tuple(pi[j][x[sp[j]]] for j in range(N))] for x in X])
            assert np.allclose(MOD[img], MOD)
            perms.append(img)
    return perms


def orbits(signs, perms):
    """Orbits of Aut(S60) on the real points (lists of indices into signs)."""
    index = {key(x): i for i, x in enumerate(signs)}
    parent = list(range(len(signs)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for i, x in enumerate(signs):
        for img in perms:
            y = np.zeros(len(X), dtype=np.int64); y[img] = x
            a, b = find(i), find(index[key(y)])
            if a != b:
                parent[a] = b
    orb = collections.defaultdict(list)
    for i in range(len(signs)):
        orb[find(i)].append(i)
    return sorted(orb.values(), key=lambda m: (len(m), min(m)))


def to_face_signs(w, tol=1e-9):
    """If w is (up to a party permutation and level relabelling) a real state supported on S60 with the face moduli,
    return its sign pattern on S60; else None."""
    supp = np.flatnonzero(np.abs(w) > tol)
    if len(supp) != len(X) or np.max(np.abs(w.imag if np.iscomplexobj(w) else 0)) > tol:
        return None
    w = np.real(w)
    Y = [tuple(int(c) for c in np.unravel_index(i, (D,) * N)) for i in supp]
    amp = dict(zip(Y, w[supp]))
    for sp in itertools.permutations(range(N)):
        sols = relabelings(Y, set(X), sp)
        if sols:
            pi = sols[0]
            vec = np.zeros(len(X))
            for y in Y:
                vec[POS[tuple(pi[j][y[sp[j]]] for j in range(N))]] = amp[y]
            if np.allclose(vec ** 2 * 72, MOD):
                return (vec < 0).astype(np.int64)
    return None
