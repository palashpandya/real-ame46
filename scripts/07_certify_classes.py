"""Rigorous separation of the 20 real classes (needs python-flint; interval arithmetic with Arb).

Invariant: N(t) = number of eigenvalues < t of 72 J^T J, J the Jacobian of psi -> (rho_A) over all six 2-party
marginals.  Local unitaries act orthogonally on the domain and on the target, party permutations permute the
marginals, so the spectrum -- hence N(t) -- is invariant under both.

Certificate.  For a real state with sqrt72 psi = a + b sqrt2 (a, b integer vectors), 72 J^T J = P + sqrt2 R with
integer P, R (exact).  It splits into blocks along the connected components of its sparsity graph.  For each block
G and threshold t, take floating-point eigenvectors V and form X = V^T (G - t I) V in Arb ball arithmetic.  If every
Gershgorin interval of X excludes 0, then X (hence V) is nonsingular and, by Sylvester's law of inertia, the number of
negative intervals equals the number of eigenvalues of G below t.  This is a proof, modulo the correctness of Arb.

Result: the pair (N(149/8), N(3/2)) takes 20 distinct values on the 20 class representatives of both faces, so the
classes are pairwise LU-inequivalent (also up to party permutations).  N(1/2) = 141 + first-order deformation dim.
"""
from fractions import Fraction

import numpy as np
from scipy.sparse.csgraph import connected_components

from _common import check, done
from ame46real import core
from ame46real.jacobian import _jac_block
from ame46real.realpoints import Face

try:
    from flint import arb, arb_mat, ctx
except ImportError:
    raise SystemExit("python-flint is required:  pip install python-flint")
ctx.prec = 200
SQRT2 = arb(2).sqrt()
THRESHOLDS = (Fraction(149, 8), Fraction(3, 2), Fraction(1, 2))


def exact_gram(F, x):
    """72 J^T J = P + sqrt2 R for the real state of face F with sign pattern x (sqrt72 psi = a + b sqrt2)."""
    a, b = np.zeros(core.DIM), np.zeros(core.DIM)
    sgn = 1 - 2 * np.asarray(x)
    for i, flat in enumerate(F.FLAT):
        (b if F.officer[i] else a)[flat] = sgn[i]
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


classes = core.load_json("real_classes.json")["classes"]
faces = {s: Face(s) for s in "+-"}
counts = {}
print("class face  N(149/8)  N(3/2)  N(1/2)")
for c in classes:
    F = faces[c["face"]]
    P, R = exact_gram(F, F.flips(c["flip"]))
    counts[c["cls"]] = tuple(certified_count(P, R, t) for t in THRESHOLDS)
    print(f"{c['cls']:>5}  {c['face']}    {counts[c['cls']][0]!s:>8}  {counts[c['cls']][1]!s:>6}  {counts[c['cls']][2]!s:>6}",
          flush=True)
check("every count certified (no Gershgorin interval touches 0)", all(None not in v for v in counts.values()))
check("certified counts equal the class table", all(counts[c["cls"]][:2] == (c["N_149_8"], c["N_3_2"]) for c in classes))
check("(N(149/8), N(3/2)) takes 20 distinct values -> 20 pairwise LU-inequivalent classes",
      len({v[:2] for v in counts.values()}) == 20)
check("N(1/2) = 141 + deformation dim (155 on the gluing locus, 152 off it)",
      all(counts[c["cls"]][2] == 141 + c["deformation_dim"] for c in classes))
done()
