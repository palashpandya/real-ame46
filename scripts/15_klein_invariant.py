"""The degree-8 Klein invariant K of Hontarenko and Zyczkowski on the 20 real classes.

For T = 6 psi (indices a = (a0, a1, a2, a3) etc.), Hontarenko and Zyczkowski (unpublished, announced October 2026) define
    K(T) = sum_{a,b,c,d} T_a T_b T_c T_d  conj(T[a0 b1 c2 d3]) conj(T[b0 a1 d2 c3]) conj(T[c0 d1 a2 b3]) conj(T[d0 c1 b2 a3]).
For every party the four conjugated factors take the indices of a, b, c, d in a permutation (identity, (ab)(cd),
(ac)(bd), (ad)(bc): the Klein four-group), so each local unitary cancels against its conjugate and K is invariant
under U(6)^4; relabelling the dummy indices makes it invariant under permutations of the parties as well.
Checks:
  (1) the contraction is invariant under Haar-random local unitaries and all 24 party permutations (numerical);
  (2) exact values in Z[sqrt2] on the 20 class representatives: on the 2T scale every non-zero coefficient is +-2
      (officers) or +-sqrt2, so each of the terms that survive on the support is +-2^m or +-2^m sqrt2, and
      K(T) = K(2T) / 2^8.  K = 72 on classes 1-2, 40 on classes 3-8 and K <= 0 on classes 9-20;
  (3) the exact values agree with a dense floating-point contraction;
  (4) K takes ten distinct values on the twenty classes, so it does not separate them (the certified invariants of
      07 do); the value 108 reported for the 96-term real state of Hontarenko and Zyczkowski is not among them, so
      that state is not equivalent to any of the 4096 real states under local unitaries and party permutations.
"""
import itertools

import numpy as np

from _common import check, done
from ame46real import core
from ame46real.realpoints import Face

D = core.D
REPORTED_96 = 108                      # K of the 96-term real state of Hontarenko and Zyczkowski (their Table 1)


def klein_dense(T):
    """K(T) by pairwise contraction (cost 6^10)."""
    Tc = np.conj(T)
    R12 = np.einsum("BCDFKP,FGHBOL->CDKPGHOL", np.einsum("ABCD,AFKP->BCDFKP", T, Tc),
                    np.einsum("EFGH,EBOL->FGHBOL", T, Tc))
    R34 = np.einsum("JKLNCH,NOPJGD->KLCHOPGD", np.einsum("IJKL,INCH->JKLNCH", T, Tc),
                    np.einsum("MNOP,MJGD->NOPJGD", T, Tc))
    return np.einsum("CDKPGHOL,KLCHOPGD->", R12, R34)


def klein_terms(flat):
    """Index octuples (a, b, c, d, X1, X2, X3, X4) into the support list `flat` of the terms of K that survive on it."""
    flat = np.asarray(flat)
    pos = -np.ones(D ** 4, dtype=np.int64)
    pos[flat] = np.arange(len(flat))
    dig = np.array(np.unravel_index(flat, (D,) * 4)).T
    n = len(flat)
    b, c, d = (x.ravel() for x in np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij"))
    out = []
    for a in range(n):
        A, B, C, Dd = dig[a], dig[b], dig[c], dig[d]
        X = [pos[np.ravel_multi_index(x, (D,) * 4)] for x in
             ((A[0], B[:, 1], C[:, 2], Dd[:, 3]), (B[:, 0], A[1], Dd[:, 2], C[:, 3]),
              (C[:, 0], Dd[:, 1], A[2], B[:, 3]), (Dd[:, 0], C[:, 1], B[:, 2], A[3]))]
        ok = np.all([x >= 0 for x in X], axis=0)
        out.append(np.stack([np.full(ok.sum(), a), b[ok], c[ok], d[ok]] + [x[ok] for x in X], axis=1))
    return np.concatenate(out)


def klein_exact(F, x, terms):
    """K(6 psi) exactly, psi the real state of face F with sign pattern x; returns (integer part, sqrt2 part)."""
    sgn = 1 - 2 * np.asarray(x)[terms]                       # signs of the 8 factors of each term
    off = F.officer[terms]                                   # factor = +-2 (officer) or +-sqrt2 on the 2T scale
    n_off = off.sum(axis=1)
    n_rt2 = 8 - n_off
    power = n_off + n_rt2 // 2                               # |term| = 2^power, times sqrt2 if n_rt2 is odd
    s = np.prod(sgn, axis=1)
    whole = sum(int(v) << int(p) for v, p in zip(s[n_rt2 % 2 == 0], power[n_rt2 % 2 == 0]))
    root2 = sum(int(v) << int(p) for v, p in zip(s[n_rt2 % 2 == 1], power[n_rt2 % 2 == 1]))
    assert whole % 2 ** 8 == 0 and root2 % 2 ** 8 == 0
    return whole // 2 ** 8, root2 // 2 ** 8


rng = np.random.default_rng(15)
table = core.load_json("real_classes.json")["classes"]
faces = {s: Face(s) for s in "+-"}


def haar(n):
    q, r = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
    return q * (np.diag(r) / np.abs(np.diag(r)))


# (1) invariance of the contraction
dev = 0.0
for r in (table[0], table[12]):
    F = faces[r["face"]]
    T = 6 * F.state(F.flips(r["flip"])).reshape((D,) * 4)
    K0 = klein_dense(T)
    for perm in itertools.permutations(range(4)):
        dev = max(dev, abs(klein_dense(np.transpose(T, perm)) - K0))
    for _ in range(3):
        U = [haar(D) for _ in range(4)]
        TU = np.einsum("ai,bj,ck,dl,ijkl->abcd", *U, T)
        dev = max(dev, abs(klein_dense(np.transpose(TU, rng.permutation(4))) - K0))
check("K is invariant under all 24 party permutations and Haar-random local unitaries (numerical)", dev < 1e-9,
      f"max deviation {dev:.1e}")

# (2)-(3) exact values on the 20 representatives
terms = {s: klein_terms(F.FLAT) for s, F in faces.items()}
values, ok_exact, ok_float = {}, True, 0.0
for r in table:
    F = faces[r["face"]]
    x = F.flips(r["flip"])
    whole, root2 = klein_exact(F, x, terms[r["face"]])
    ok_exact &= root2 == 0
    values[r["cls"]] = whole
    ok_float = max(ok_float, abs(klein_dense(6 * F.state(x).reshape((D,) * 4)) - whole))
check("K(6 psi) is an integer on each of the 20 representatives (exact, Z[sqrt2])", ok_exact,
      f"{len(terms['+'])} and {len(terms['-'])} surviving terms on the two supports")
check("exact values agree with the dense floating-point contraction", ok_float < 1e-9, f"max deviation {ok_float:.1e}")
check("K = 72 on classes 1-2, 40 on classes 3-8, K <= 0 on classes 9-20",
      all(values[c] == 72 for c in (1, 2)) and all(values[c] == 40 for c in range(3, 9))
      and all(values[c] <= 0 for c in range(9, 21)))

# (4) separation
distinct = sorted(set(values.values()), reverse=True)
check("K takes 10 distinct values on the 20 classes", len(distinct) == 10, f"{distinct}")
check(f"K = {REPORTED_96} (the 96-term real state of Hontarenko and Zyczkowski) is not among them",
      REPORTED_96 not in distinct)

print("\nclass face  flipped        K(6 psi)")
for r in table:
    print(f"{r['cls']:>5}  {r['face']}    {' '.join(r['flip']) or '-':<14} {values[r['cls']]:>8}")
done()
