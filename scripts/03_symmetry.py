"""The order-4 local symmetry behind the family.

The state was found by searching for AME states inside the -1 eigenspace of a local symmetry U^{(x)4} with
U = V diag(i, -1, -1, -i, -i, -i) V^dag (V Haar random; the spectrum was selected by necessary conditions, using only
N = 4 and d = 6).  In the frame of this repository the symmetry is diagonal:  u = diag(i, -1, -1, -i, -i, -i) on every
party, levels 0..5.  Checks (exact unless stated):
  (1) the -1 eigenspace of u^{(x)4} is spanned by the 360 basis states x with e(x_1)+...+e(x_4) = 2 (mod 4),
      e = (1, 2, 2, 3, 3, 3);
  (2) all 72 points of L u O u N lie in it, so EVERY state of both families -- in particular the real state psi and
      the uniform-amplitude states -- satisfies u^{(x)4}|phi> = -|phi>;
  (3) the diagonal local unitaries acting as a scalar on the support form (mod trivial ones) the group listed below;
      u^{(x)4} has order 4 in it;
  (4) the diagonal symmetries of this spectral type (every party {i,-1,-1,-i,-i,-i}, eigenvalue -1) of the support:
      u itself and its images under relabelling the levels {0,3,4,5}.
The commutant of u is U(1) x U(2) x U(3) (levels {0}, {1,2}, {3,4,5}); frame changes inside it keep the sector, and
minimising the support over them (and over U(6)) turned the search hit (support 348) into the support-72 form.
"""
import itertools

import numpy as np

from _common import check, done
from ame46real import core, family as fam
from ame46real.lattice import invariant_factors

e = np.array([1, 2, 2, 3, 3, 3])                                # u = diag(i^e)
sector = [x for x in itertools.product(range(6), repeat=4) if e[list(x)].sum() % 4 == 2]
check("(1) dim of the -1 eigenspace of u^(x)4", len(sector) == 360, f"{len(sector)}")
S72 = np.array(fam.POINTS)
check("(2) all 72 support points lie in the sector", all(e[x].sum() % 4 == 2 for x in S72))
u = np.diag(1j ** e)
U4 = lambda v: np.einsum("ai,bj,ck,dl,ijkl->abcd", u, u, u, u, v.reshape((6,) * 4)).ravel()
rng = np.random.default_rng(4)
tests = {"psi (real)": core.psi(),
         "psi_sigma(delta, theta), random": fam.psi72(1, rng.uniform(-1, 1), rng.uniform(0, 7, 9)),
         "psi_face(theta), random": fam.psi_face(rng.uniform(0, 7, 11)),
         "uniform psi_-(0, 0)": fam.psi72(-1, 0.0, np.zeros(9))}
for name, v in tests.items():
    check(f"(2) u^(x)4 |phi> = -|phi> for {name} (numerical)", np.allclose(U4(v), -v))


def scalar_group(qs):
    """Invariant factors > 1 of [G | -1]: the finite group of diagonal local unitaries that act on every point of the
    support by the same phase, modulo the trivial ones (independent phases per party)."""
    G = fam.gauge_matrix(qs)
    B = np.hstack([G, -np.ones((len(qs), 1), dtype=np.int64)])
    assert np.linalg.matrix_rank(B) == 21                       # continuous part = per-party phases only
    return [f for f in invariant_factors(B.tolist()) if f > 1]


g72, g60 = scalar_group(fam.S72), scalar_group(fam.S60)
check("(3) diagonal symmetry group of the support-72 array", sorted(g72) == [2, 2, 4], " x ".join(f"Z{f}" for f in g72))
check("(3) diagonal symmetry group of the support-60 face", True, " x ".join(f"Z{f}" for f in g60))
check("(3) u^(x)4 has order 4 in this group (u, u^2, u^3 are not scalar on a party; u^4 = 1)",
      all(len(set((k * e) % 4)) > 1 for k in (1, 2, 3)) and not np.any((4 * e) % 4))
sols = []
for a in itertools.product(*[sorted(set(itertools.permutations((1, 2, 2, 3, 3, 3))))] * 3):
    a = [np.array(t) for t in a]
    need = (2 - a[0][S72[:, 0]] - a[1][S72[:, 1]] - a[2][S72[:, 2]]) % 4
    e3 = {}
    if all(e3.setdefault(v, n) == n for v, n in zip(S72[:, 3], need)) and len(e3) == 6 \
            and sorted(e3.values()) == [1, 2, 2, 3, 3, 3]:
        sols.append(tuple(tuple(int(c) for c in t) for t in a) + (tuple(int(e3[v]) for v in range(6)),))
check("(4) diagonal symmetries of this spectral type containing the support: 4, all with -1 on levels {1,2};"
      " they differ only in which level of {0,3,4,5} carries i (u is one of them)",
      len(sols) == 4 and tuple(tuple(e.tolist()) for _ in range(4)) in sols
      and all(t[1] == t[2] == 2 for sol in sols for t in sol))
done()
