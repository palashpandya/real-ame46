"""Real AME states on the two support-60 faces of the support-72 family, and their monomial symmetries.

Face "+" : support L u O, officers O (coefficient 1/6), base state psi       (delta = +1 end of the family).
Face "-" : support L u N, officers N (coefficient 1/6), base state psi_minus (delta = -1 end of the family).
A real state with support exactly the face has signs x in GF(2)^60 (1 = minus) with Mg x = 1, where Mg is the
block-congruence matrix mod 2 (officers are 1x1 blocks, so only the light 2x2 blocks give congruences).  Two sign
patterns are the same state up to diagonal local phases iff Pc (x - x') = 0 (mod 2), where Pc are the rows of the
Smith transform P of the gauge matrix beyond its rank.
"""
import collections
import itertools

import numpy as np

from .core import D, N, DIM
from . import family as fam
from .lattice import snf, gf2_nullspace, gf2_rank


def _cokernel_rows(B):
    D_, P, _ = snf(np.asarray(B).tolist())
    r = sum(1 for i in range(min(len(D_), len(D_[0]))) if D_[i][i])
    return np.array(P, dtype=np.int64)[r:]


def relabelings(Y, target, sp, first_only=False):
    """Level relabellings pi (one per party) with {pi(Y[:, sp])} = target, found by backtracking."""
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
        if first_only and out:
            return
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


class Face:
    def __init__(self, sign):
        assert sign in "+-"
        self.sign = sign
        off_type = "O" if sign == "+" else "N"
        self.qs = [q for q in fam.S72 if fam.TYPE[q] in ("L", off_type)]
        self.X = [fam.POINTS[q] for q in self.qs]
        self.n = len(self.X)
        self.POS = {x: i for i, x in enumerate(self.X)}
        self.FLAT = [fam.FLAT[q] for q in self.qs]
        self.officer = np.array([fam.TYPE[q] == off_type for q in self.qs])
        self.MOD = np.where(self.officer, 2.0, 1.0)                        # 72 |amplitude|^2
        self.X0 = np.array([int(fam.EPS[q] < 0) if fam.TYPE[q] == "L" else 0 for q in self.qs])
        self.Mg = fam.congruence_matrix(self.qs)[0] % 2
        G = fam.gauge_matrix(self.qs)
        self.Prest = _cokernel_rows(G)                                    # integer rows
        self.Pc = self.Prest % 2
        self.P2c = _cokernel_rows(np.hstack([fam.K[self.qs], G])) % 2      # (x - x0) on the gluing locus

    def state(self, x):
        v = np.zeros(DIM)
        v[self.FLAT] = np.sqrt(self.MOD / 72) * (1 - 2 * np.asarray(x))
        return v

    def key(self, x):
        return tuple((self.Pc @ np.asarray(x)) % 2)

    def flips(self, labels):
        """Sign pattern of the base state with the given officer coefficients negated."""
        x = self.X0.copy()
        for s in labels:
            i = self.POS[tuple(int(c) for c in s)]
            assert self.officer[i]
            x[i] ^= 1
        return x

    def real_points(self):
        """(sign patterns of the distinct real points, dim of the GF(2) solution space, dim of its gauge part)."""
        assert not np.any((self.Mg @ self.X0 - 1) % 2)
        ker = gf2_nullspace(self.Mg)
        comp, img = [], []
        for v in ker:
            trial = img + [(self.Pc @ v) % 2]
            if gf2_rank(np.array(trial)) == len(trial):
                comp.append(v); img = trial
        signs = [(self.X0 + sum(c * v for c, v in zip(cs, comp))) % 2
                 for cs in itertools.product((0, 1), repeat=len(comp))]
        return signs, len(ker), len(ker) - len(comp)

    def on_C(self, x):
        """Is the real state a limit of support-72 AME states (the base state being one)?"""
        return not np.any((self.P2c @ ((np.asarray(x) - self.X0) % 2)) % 2)

    def real_up_to_phases(self, turns):
        """Exact: are the phases (turns, on this face) those of a real state up to diagonal local unitaries?"""
        from fractions import Fraction
        y = [sum(int(p) * Fraction(t) for p, t in zip(row, turns)) for row in self.Prest]
        return all((2 * v).denominator == 1 for v in y)

    def automorphisms(self):
        """Party permutations x level relabellings mapping the support onto itself and officers onto officers,
        as permutations of the 60 points (img[i] = image of point i)."""
        S = set(self.X)
        perms = []
        for sp in itertools.permutations(range(N)):
            for pi in relabelings(self.X, S, sp):
                img = np.array([self.POS[tuple(pi[j][x[sp[j]]] for j in range(N))] for x in self.X])
                if np.array_equal(self.officer[img], self.officer):
                    perms.append(img)
        return perms

    def orbits(self, signs, perms):
        """Orbits of the automorphism group on the real points (lists of indices into signs)."""
        index = {self.key(x): i for i, x in enumerate(signs)}
        parent = list(range(len(signs)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]; a = parent[a]
            return a
        for i, x in enumerate(signs):
            for img in perms:
                y = np.zeros(self.n, dtype=np.int64); y[img] = x
                a, b = find(i), find(index[self.key(y)])
                if a != b:
                    parent[a] = b
        orb = collections.defaultdict(list)
        for i in range(len(signs)):
            orb[find(i)].append(i)
        return sorted(orb.values(), key=lambda m: (len(m), min(m)))

    def to_face_signs(self, w, tol=1e-9):
        """If the real vector w is, up to a party permutation and a level relabelling, a state supported on this
        face with officers on officers and the face moduli, return its sign pattern; else None."""
        w = np.asarray(w)
        supp = np.flatnonzero(np.abs(w) > tol)
        if len(supp) != self.n or (np.iscomplexobj(w) and np.max(np.abs(w.imag)) > tol):
            return None
        w = np.real(w)
        Y = [tuple(int(c) for c in np.unravel_index(i, (D,) * N)) for i in supp]
        amp = dict(zip(Y, w[supp]))
        for sp in itertools.permutations(range(N)):
            for pi in relabelings(Y, set(self.X), sp):
                vec = np.zeros(self.n)
                for y in Y:
                    vec[self.POS[tuple(pi[j][y[sp[j]]] for j in range(N))]] = amp[y]
                if np.allclose(vec ** 2 * 72, self.MOD):
                    return (vec < 0).astype(np.int64)
        return None


def monomially_equivalent(f1, f2):
    """Is there a party permutation x level relabelling mapping support f1 onto support f2, officers onto officers?"""
    for sp in itertools.permutations(range(N)):
        for pi in relabelings(f1.X, set(f2.X), sp):
            img = [f2.POS[tuple(pi[j][x[sp[j]]] for j in range(N))] for x in f1.X]
            if np.array_equal(f2.officer[img], f1.officer):
                return True
    return False


# module-level API for face "+"
_PLUS = Face("+")
X, POS, FLAT, MOD, X0 = _PLUS.X, _PLUS.POS, _PLUS.FLAT, _PLUS.MOD, _PLUS.X0
state, key, on_C, to_face_signs = _PLUS.state, _PLUS.key, _PLUS.on_C, _PLUS.to_face_signs


def real_points():
    return _PLUS.real_points()


def automorphisms():
    return _PLUS.automorphisms()


def orbits(signs, perms):
    return _PLUS.orbits(signs, perms)


def theta_points():
    """Distinct real points of face "+" reached by theta_j in {0, pi} in the face closed form."""
    Kf = fam.face_K()
    return {key((X0 + Kf @ np.array(t)) % 2) for t in itertools.product((0, 1), repeat=11)}
