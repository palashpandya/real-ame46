"""Exact integer linear algebra: Smith normal form, solving A y = r (mod 1) over Q, and GF(2) linear algebra."""
from fractions import Fraction

import numpy as np


def snf(A):
    """Smith normal form with transforms: returns D, P, Q (lists of Python ints) with P A Q = D."""
    A = [[int(v) for v in row] for row in A]
    m, n = len(A), len(A[0])
    P = [[int(i == j) for j in range(m)] for i in range(m)]
    Q = [[int(i == j) for j in range(n)] for i in range(n)]

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]; P[i], P[j] = P[j], P[i]

    def swap_cols(i, j):
        for R in A: R[i], R[j] = R[j], R[i]
        for R in Q: R[i], R[j] = R[j], R[i]

    def add_row(dst, src, k):
        A[dst] = [x + k * y for x, y in zip(A[dst], A[src])]; P[dst] = [x + k * y for x, y in zip(P[dst], P[src])]

    def add_col(dst, src, k):
        for R in A: R[dst] += k * R[src]
        for R in Q: R[dst] += k * R[src]

    t = 0
    while t < min(m, n):
        nz = [(abs(A[i][j]), i, j) for i in range(t, m) for j in range(t, n) if A[i][j]]
        if not nz:
            break
        _, i, j = min(nz); swap_rows(t, i); swap_cols(t, j)
        while True:
            changed = False
            for i in range(t + 1, m):
                if A[i][t]:
                    add_row(i, t, -(A[i][t] // A[t][t]))
                    changed |= bool(A[i][t])
            for j in range(t + 1, n):
                if A[t][j]:
                    add_col(j, t, -(A[t][j] // A[t][t]))
                    changed |= bool(A[t][j])
            if changed:
                nz = [(abs(A[i][j]), i, j) for i in range(t, m) for j in range(t, n)
                      if A[i][j] and (i == t or j == t)]
                _, i, j = min(nz); swap_rows(t, i); swap_cols(t, j)
                continue
            bad = next(((i, j) for i in range(t + 1, m) for j in range(t + 1, n) if A[i][j] % A[t][t]), None)
            if bad is None:
                break
            add_row(t, bad[0], 1)
        if A[t][t] < 0:
            A[t] = [-x for x in A[t]]; P[t] = [-x for x in P[t]]
        t += 1
    return A, P, Q


def invariant_factors(A):
    D, _, _ = snf(A)
    return [D[i][i] for i in range(min(len(D), len(D[0]))) if D[i][i]]


def solve_mod1(A, r):
    """A rational solution y of A y = r (mod 1), or None.  A: integer matrix, r: rationals (turns)."""
    D, P, Q = snf(A)
    m, n = len(A), len(A[0])
    rk = sum(1 for i in range(min(m, n)) if D[i][i])
    Pr = [sum(Fraction(P[i][j]) * Fraction(r[j]) for j in range(m)) for i in range(m)]
    if any(Pr[i].denominator != 1 for i in range(rk, m)):
        return None
    z = [Pr[i] / D[i][i] if i < rk else Fraction(0) for i in range(n)]
    return [sum(Fraction(Q[i][j]) * z[j] for j in range(n)) for i in range(n)]


def all_solutions_mod1(A, r):
    """All solutions y in (Q/Z)^n of A y = r (mod 1) when the solution set is finite (A has full column rank)."""
    import itertools
    D, P, Q = snf(A)
    m, n = len(A), len(A[0])
    inv = [D[i][i] for i in range(min(m, n)) if D[i][i]]
    assert len(inv) == n, "solution set is not finite"
    Pr = [sum(Fraction(P[i][j]) * Fraction(r[j]) for j in range(m)) for i in range(m)]
    if any(Pr[i].denominator != 1 for i in range(n, m)):
        return []
    sols = []
    for ks in itertools.product(*[range(abs(v)) for v in inv]):
        z = [(Pr[i] + ks[i]) / inv[i] for i in range(n)]
        sols.append([sum(Fraction(Q[i][j]) * z[j] for j in range(n)) % 1 for i in range(n)])
    return sols


def gf2_rref(M):
    A = np.array(M, dtype=np.int64) % 2
    r, piv = 0, []
    for c in range(A.shape[1]):
        p = next((i for i in range(r, A.shape[0]) if A[i, c]), None)
        if p is None:
            continue
        A[[r, p]] = A[[p, r]]
        for i in range(A.shape[0]):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        piv.append(c); r += 1
    return A[:r], piv


def gf2_rank(M):
    return len(gf2_rref(M)[1])


def gf2_nullspace(M):
    M = np.array(M, dtype=np.int64)
    B, piv = gf2_rref(M)
    out = []
    for f in [c for c in range(M.shape[1]) if c not in piv]:
        v = np.zeros(M.shape[1], dtype=np.int64); v[f] = 1
        for row, p in zip(B, piv):
            if row[f]:
                v[p] = 1
        out.append(v)
    return np.array(out, dtype=np.int64).reshape(len(out), M.shape[1])


def gf2_solvable(M, b):
    """Is M x = b solvable over GF(2)?"""
    return gf2_rank(np.column_stack([M, b])) == gf2_rank(M)
