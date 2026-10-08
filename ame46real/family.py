"""The support-72 family psi_sigma(delta, theta) and the support-60 face psi_face(theta) (closed forms).

Points are basis states of the 72-point support L u O u N (data/family.json), in lexicographic order.
Phases are handled exactly in *turns* (multiples of 2 pi) with Fractions.
"""
from fractions import Fraction

import numpy as np

from .core import DIM, D, N, index, load_json, tup, two_by_two_blocks
from .lattice import snf, solve_mod1

F = load_json("family.json")
POINTS = [tup(p["x"]) for p in F["points"]]
TYPE = [p["type"] for p in F["points"]]
EPS = np.array([p["eps"] for p in F["points"]])
K = np.array([p["K"] for p in F["points"]], dtype=np.int64)            # 72 x 9 phase directions
FLAT = [index(x) for x in POINTS]
POS = {x: q for q, x in enumerate(POINTS)}
S72 = list(range(len(POINTS)))
S60 = [q for q in S72 if TYPE[q] != "N"]                                # support of psi (L u O)
NPTS = [q for q in S72 if TYPE[q] == "N"]


def face_K():
    """60 x 11 phase directions on the face: K restricted to L u O, plus theta_10 and theta_11 (single officers)."""
    Kf = np.zeros((len(S60), 11), dtype=np.int64)
    Kf[:, :9] = K[S60]
    for j, s in enumerate(F["face_extra"]):
        Kf[S60.index(POS[tup(s)]), 9 + j] = 1
    return Kf


def psi72(sigma, delta, theta):
    """Eq. (F72): (1/sqrt72)[sum_L eps e^{i K theta} + sqrt(1+delta) sum_O e^{i K theta} + i sigma sqrt(1-delta) sum_N e^{i K theta}]."""
    amp = np.array([EPS[q] if t == "L" else (np.sqrt(1 + delta) if t == "O" else 1j * sigma * np.sqrt(1 - delta))
                    for q, t in enumerate(TYPE)], dtype=complex)
    v = np.zeros(DIM, dtype=complex)
    v[FLAT] = amp * np.exp(1j * (K @ np.asarray(theta, float))) / np.sqrt(72)
    return v


def psi_face(theta):
    """Eq. (F60): (1/sqrt72) sum_L eps e^{i K' theta} + (1/6) sum_O e^{i K' theta}."""
    amp = np.array([EPS[q] * (1.0 if TYPE[q] == "L" else np.sqrt(2)) for q in S60], dtype=complex)
    v = np.zeros(DIM, dtype=complex)
    v[[FLAT[q] for q in S60]] = amp * np.exp(1j * (face_K() @ np.asarray(theta, float))) / np.sqrt(72)
    return v


# ---------- exact phase bookkeeping (turns) ------------------------------------------------------------------
def base_turns(sigma):
    """Phases (turns) of psi_sigma(delta, 0) on the 72 points."""
    return [Fraction(0) if (t != "N" and EPS[q] > 0) else
            (Fraction(1, 2) if t != "N" else (Fraction(1, 4) if sigma > 0 else Fraction(3, 4)))
            for q, t in enumerate(TYPE)]


def turns72(sigma, theta_over_pi):
    """Exact phases (turns) of psi_sigma(delta, pi * t) for rational t."""
    t = [Fraction(v) for v in theta_over_pi]
    return [(b + sum(int(K[q, j]) * t[j] for j in range(9)) / 2) % 1 for q, b in enumerate(base_turns(sigma))]


def congruence_matrix(qs):
    """Rows phi_a - phi_b - phi_c + phi_e for every 2x2 block [[a, b], [c, e]] of the support {POINTS[q] : q in qs}.
    AME on such a support (with equal moduli on block diagonals) <=> M phi = 1/2 (mod 1) in turns."""
    pts = [POINTS[q] for q in qs]
    pos = {x: i for i, x in enumerate(pts)}
    quads, singles = two_by_two_blocks(pts)
    M = np.zeros((len(quads), len(qs)), dtype=np.int64)
    for r, (a, b, c, e) in enumerate(quads):
        M[r, pos[a]] += 1; M[r, pos[b]] -= 1; M[r, pos[c]] -= 1; M[r, pos[e]] += 1
    return M, quads, singles


def gauge_matrix(qs):
    """Diagonal local phases: G[i, 6 j + v] = 1 if point i has level v at party j."""
    return np.array([[int(POINTS[q][j] == v) for j in range(N) for v in range(D)] for q in qs], dtype=np.int64)


def equal_up_to_local_phases(t1, t2, qs):
    """Do two phase patterns (turns, on the points qs) differ by diagonal local unitaries?  Exact."""
    diff = [(Fraction(a) - Fraction(b)) % 1 for a, b in zip(t1, t2)]
    return solve_mod1(gauge_matrix(qs).tolist(), diff) is not None


class Components:
    """Exact component test for the solution torus {t : M t = 1/2 (mod 1)}: t1, t2 lie in the same connected
    component iff t1 - t2 is in ker_R(M) + Z^n, i.e. iff M (t1 - t2) lies in the lattice M Z^n."""

    def __init__(self, M):
        self.M = np.asarray(M, dtype=np.int64)
        D_, self.P, _ = snf(self.M.tolist())
        self.d = [D_[i][i] for i in range(min(len(D_), len(D_[0]))) if D_[i][i]]
        self.count = int(np.prod(self.d))

    def same(self, t1, t2):
        dl = [Fraction(a) - Fraction(b) for a, b in zip(t1, t2)]
        Md = [sum(int(self.M[r, i]) * dl[i] for i in range(len(dl))) for r in range(self.M.shape[0])]
        assert all(v.denominator == 1 for v in Md), "not two solutions of the same congruences"
        w = [sum(self.P[i][j] * int(Md[j]) for j in range(len(Md))) for i in range(len(self.P))]
        return all(w[i] % self.d[i] == 0 for i in range(len(self.d)))
