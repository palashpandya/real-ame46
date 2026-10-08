"""Basic objects: the four-quhex Hilbert space, cuts, the state psi and the AME test."""
import itertools
import json
import os

import numpy as np

D, N = 6, 4
DIM = D ** N
CUTS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]       # the three 2|2 bipartitions
PAIRS = list(itertools.combinations(range(N), 2))                      # the six 2-party marginals
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "data")

# coefficients of sqrt(72) psi as a + b sqrt2
COEFF = {"+1": (1, 0), "-1": (-1, 0), "+sqrt2": (0, 1), "-sqrt2": (0, -1)}


def load_json(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


def tup(label):
    return tuple(int(c) for c in label)


def lab(x):
    return "".join(str(int(c)) for c in x)


def index(x):
    return int(np.ravel_multi_index(tuple(x), (D,) * N))


def psi_terms():
    """{x: (a, b)} with sqrt(72) psi_x = a + b sqrt2."""
    return {tup(x): COEFF[c] for x, c in load_json("psi_terms.json")["terms"]}


def psi():
    v = np.zeros(DIM)
    for x, (a, b) in psi_terms().items():
        v[index(x)] = (a + b * np.sqrt(2)) / np.sqrt(72)
    return v


def flip(v, labels):
    """Negate the coefficients of the given basis states."""
    w = np.array(v, copy=True)
    for s in labels:
        w[index(tup(s))] *= -1
    return w


def cut_matrix(v, cut):
    A, B = cut
    return np.transpose(np.asarray(v).reshape((D,) * N), A + B).reshape(D * D, D * D)


def ame_error(v):
    """max over the three cuts of max |singular value - 1/6|  (0 for an AME(4,6) state)."""
    v = np.asarray(v) / np.linalg.norm(v)
    return max(float(np.max(np.abs(np.linalg.svd(cut_matrix(v, c), compute_uv=False) - 1 / D))) for c in CUTS)


def support(v, tol=1e-9):
    return [tuple(int(c) for c in np.unravel_index(i, (D,) * N)) for i in np.flatnonzero(np.abs(v) > tol)]


def point(row, col, cut):
    """Basis state with row label `row` (pair of levels) and column label `col` in the given cut."""
    A, B = cut
    x = [0] * N
    x[A[0]], x[A[1]], x[B[0]], x[B[1]] = row[0], row[1], col[0], col[1]
    return tuple(x)


def blocks(points, cut):
    """Connected components of the row/column graph of a support in one cut.
    Returns a list of (rows, cols) with sorted row and column labels (pairs of levels)."""
    A, B = cut
    pts = set(points)
    rows = {}
    for x in pts:
        rows.setdefault((x[A[0]], x[A[1]]), set()).add((x[B[0]], x[B[1]]))
    cols = {}
    for r, cs in rows.items():
        for c in cs:
            cols.setdefault(c, set()).add(r)
    seen, out = set(), []
    for r0 in sorted(rows):
        if r0 in seen:
            continue
        R, C, stack = {r0}, set(), [("r", r0)]
        while stack:
            kind, v = stack.pop()
            for w in (rows[v] if kind == "r" else cols[v]):
                S = C if kind == "r" else R
                if w not in S:
                    S.add(w)
                    stack.append(("c" if kind == "r" else "r", w))
        seen |= R
        out.append((sorted(R), sorted(C)))
    return out


def two_by_two_blocks(points):
    """All 2x2 blocks [[a, b], [c, e]] (as point quadruples) of a support in the three cuts, plus the 1x1 blocks.
    Raises if a block of another shape occurs."""
    quads, singles = [], []
    for cut in CUTS:
        for R, C in blocks(points, cut):
            if len(R) == len(C) == 1:
                singles.append(point(R[0], C[0], cut))
            elif len(R) == len(C) == 2:
                quads.append(tuple(point(r, c, cut) for r in R for c in C))      # a, b, c, e
            else:
                raise ValueError(f"block of shape {len(R)}x{len(C)}")
    return quads, singles
