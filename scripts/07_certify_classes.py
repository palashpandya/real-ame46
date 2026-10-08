"""Rigorous separation of the 8 real classes (needs python-flint; interval arithmetic with Arb).

Invariant: N(t) = number of eigenvalues < t of 72 J^T J, J the Jacobian of psi -> (rho_A) over all six 2-party
marginals.  Local unitaries act orthogonally on the domain and on the target, party permutations permute the
marginals, so the spectrum -- hence N(t) -- is invariant under both.

Certificate.  For a real state with sqrt72 psi = a + b sqrt2 (a, b integer vectors), 72 J^T J = P + sqrt2 R with
integer P, R (exact).  It splits into blocks along the connected components of its sparsity graph.  For each block
G and threshold t, take floating-point eigenvectors V and form X = V^T (G - t I) V in Arb ball arithmetic.  If every
Gershgorin interval of X excludes 0, then X (hence V) is nonsingular and, by Sylvester's law of inertia, the number of
negative intervals equals the number of eigenvalues of G below t.  This is a proof, modulo the correctness of Arb.

Result: N(41/8) takes 8 distinct values on the 8 class representatives, so the classes are pairwise
LU-inequivalent (also up to party permutations); N(1/2) = 141 + first-order deformation dimension (upper bound).
"""
import numpy as np
from scipy.sparse.csgraph import connected_components

from _common import check, done
from ame46real import core
from ame46real.jacobian import _jac_block

try:
    from flint import arb, arb_mat, ctx
except ImportError:
    raise SystemExit("python-flint is required:  pip install python-flint")
ctx.prec = 200
SQRT2 = arb(2).sqrt()


def exact_gram(flipped):
    """72 J^T J = P + sqrt2 R for sqrt72 psi = a + b sqrt2 (psi with the given officer signs flipped), as integer
    sparse matrices."""
    a, b = np.zeros(core.DIM), np.zeros(core.DIM)
    for x, (p, q) in core.psi_terms().items():
        s = -1 if core.lab(x) in flipped else 1
        a[core.index(x)], b[core.index(x)] = s * p, s * q
    P = R = 0
    for A in core.PAIRS:
        Ja, Jb = _jac_block(a, A), _jac_block(b, A)
        P = P + Ja.T @ Ja + 2 * (Jb.T @ Jb)
        R = R + Ja.T @ Jb + Jb.T @ Ja
    P, R = P.tocsr(), R.tocsr()
    for Mx in (P, R):
        assert np.all(Mx.data == np.round(Mx.data))
    return P, R


def certified_count(P, R, t):
    """Number of eigenvalues of P + sqrt2 R below t (t a Fraction), certified; None if a disc touches 0."""
    n, lab = connected_components((abs(P) + abs(R)) > 0, directed=False)
    tf, ta = float(t), arb(t.numerator) / t.denominator
    total = 0
    for blk in range(n):
        ix = np.flatnonzero(lab == blk)
        Pb, Rb = P[ix][:, ix].toarray(), R[ix][:, ix].toarray()
        m = len(ix)
        _, V = np.linalg.eigh(Pb + np.sqrt(2) * Rb - tf * np.eye(m))
        G = arb_mat([[arb(int(Pb[i, j])) + arb(int(Rb[i, j])) * SQRT2 - (ta if i == j else 0) for j in range(m)]
                     for i in range(m)])
        Va = arb_mat([[arb(float(V[i, j])) for j in range(m)] for i in range(m)])
        X = Va.transpose() * G * Va
        for i in range(m):
            rad = sum((abs(X[i, j]) for j in range(m) if j != i), arb(0))
            if X[i, i] - rad > 0:
                continue
            if X[i, i] + rad < 0:
                total += 1
                continue
            return None
    return total


psi = core.psi()
classes = core.load_json("real_classes.json")["classes"]
expected = {1: 261, 2: 265, 3: 255, 4: 256, 5: 264, 6: 254, 7: 257, 8: 258}
from fractions import Fraction
counts = {}
print("class  N(41/8)  N(1/2)")
for c in classes:
    P, R = exact_gram(set(c["flip"]))
    counts[c["cls"]] = (certified_count(P, R, Fraction(41, 8)), certified_count(P, R, Fraction(1, 2)))
    print(f"{c['cls']:>5}  {counts[c['cls']][0]!s:>7}  {counts[c['cls']][1]!s:>6}", flush=True)
check("every count certified (no Gershgorin interval touches 0)", all(None not in v for v in counts.values()))
check("N(41/8) = 261, 265, 255, 256, 264, 254, 257, 258: eight distinct values -> 8 pairwise LU-inequivalent classes",
      {k: v[0] for k, v in counts.items()} == expected and len({v[0] for v in counts.values()}) == 8)
check("N(1/2) = 141 + deformation dim (155 for classes 1-2, 152 for 3-8)",
      all(counts[c["cls"]][1] == 141 + c["deformation_dim"] for c in classes))
done()
