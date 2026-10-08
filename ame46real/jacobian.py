"""First-order (tangent) analysis of AME(4,6) states in real coordinates.

Coordinates: column 2k is the real perturbation e_k, column 2k+1 the imaginary perturbation i e_k.
For a 2-party subset A, d rho_A in direction e_k is E + E^dag with E = |r><conj M[:, c]|, where M is the reshaping
A | A^c and (r, c) the row and column of k; in direction i e_k it is i (E - E^dag).
"""
import itertools

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from .core import D, DIM, N, PAIRS


def _jac_block(v, A):
    """Sparse real Jacobian of psi -> rho_A (rows: Re and Im of the 36 x 36 entries; columns as above)."""
    comp = tuple(q for q in range(N) if q not in A)
    perm = A + comp
    M = np.transpose(np.asarray(v, complex).reshape((D,) * N), perm).reshape(D * D, D * D)
    idx = np.arange(DIM).reshape((D,) * N).transpose(perm).reshape(D * D, D * D)   # idx[r, c] = basis index
    r, c, b = np.meshgrid(np.arange(D * D), np.arange(D * D), np.arange(D * D), indexing="ij")
    r, c, b = r.ravel(), c.ravel(), b.ravel()
    k = idx[r, c]
    m = M[b, c]
    # E = |r><conj M[:, c]|:  E[r, b] = conj(M[b, c]);  E^dag[b, r] = M[b, c]
    pos = np.concatenate([r * D * D + b, b * D * D + r])
    re_dir = np.concatenate([np.conj(m), m])                         # d rho for e_k:   E + E^dag
    im_dir = np.concatenate([1j * np.conj(m), -1j * m])              # d rho for i e_k: i E - i E^dag
    kk = np.concatenate([k, k])
    rows, cols, vals = [], [], []
    for colofs, z in ((0, re_dir), (1, im_dir)):
        for rowofs, part in ((0, z.real), (D ** 4, z.imag)):
            rows.append(pos + rowofs); cols.append(2 * kk + colofs); vals.append(part)
    return sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(2 * D ** 4, 2 * DIM))


def gram(v, subsets):
    """Sum over the given 2-party subsets A of J_A^T J_A (2592 x 2592, real symmetric, sparse)."""
    G = sparse.csr_matrix((2 * DIM, 2 * DIM))
    for A in subsets:
        J = _jac_block(v, A)
        G = G + J.T @ J
    return G


def eigvalsh_blocks(G):
    """Eigenvalues of a sparse symmetric matrix, computed block by block over its connected components."""
    n, lab = connected_components(abs(G) > 0, directed=False)
    out = []
    for b in range(n):
        ix = np.flatnonzero(lab == b)
        out.append(np.linalg.eigvalsh(G[ix][:, ix].toarray()))
    return np.sort(np.concatenate(out))


def lu_tangent(v):
    """Real 2592 x 144 matrix spanning the tangent space of the local-unitary orbit."""
    t = np.asarray(v, complex).reshape((D,) * N)
    cols = []
    for j in range(N):
        for a in range(D):
            for b in range(D):
                H = np.zeros((D, D), complex)
                if a == b:
                    H[a, a] = 1
                elif a < b:
                    H[a, b] = H[b, a] = 1
                else:
                    H[a, b], H[b, a] = 1j, -1j
                w = np.moveaxis(np.tensordot(1j * H, t, axes=([1], [j])), 0, j).ravel()
                cols.append(np.column_stack([w.real, w.imag]).ravel())
    return np.array(cols).T


def deformation_dim(v, tol=1e-9):
    """dim of first-order AME deformations modulo local unitaries:
    nullity of the Jacobian of (rho_A - I/36) over the three cuts, minus the rank of the LU orbit tangent.
    Returns (dim, nullity, lu_rank, gap) where gap = (largest 'zero', smallest 'nonzero') eigenvalue of J^T J."""
    v = np.asarray(v, complex) / np.linalg.norm(v)
    ev = eigvalsh_blocks(gram(v, [(0, 1), (0, 2), (0, 3)]))
    scale = ev[-1]
    nullity = int(np.sum(ev < tol * scale))
    s = np.linalg.svd(lu_tangent(v), compute_uv=False)
    rT = int(np.sum(s > 1e-9 * s[0]))
    return nullity - rT, nullity, rT, (float(ev[nullity - 1] / scale), float(ev[nullity] / scale))


def real_deformation_dim(v, tol=1e-9):
    """Real first-order deformations of a real state modulo local orthogonal maps."""
    v = np.asarray(v, float) / np.linalg.norm(v)
    G = gram(v, [(0, 1), (0, 2), (0, 3)])[0::2][:, 0::2]                 # real perturbations only
    ev = eigvalsh_blocks(G)
    nullity = int(np.sum(ev < tol * ev[-1]))
    t = v.reshape((D,) * N)
    gens = []
    for j in range(N):
        for a, b in itertools.combinations(range(D), 2):
            X = np.zeros((D, D)); X[a, b], X[b, a] = 1, -1
            gens.append(np.moveaxis(np.tensordot(X, t, axes=([1], [j])), 0, j).ravel())
    s = np.linalg.svd(np.array(gens).T, compute_uv=False)
    return nullity - int(np.sum(s > 1e-9 * s[0]))


def marginal_spectrum(v):
    """LU- and party-permutation-invariant fingerprint: eigenvalues of 72 J^T J, where J is the Jacobian of
    psi -> (rho_A)_{A in all six pairs} (= 72 x squared singular values of J)."""
    v = np.asarray(v, complex) / np.linalg.norm(v)
    return eigvalsh_blocks(72 * gram(v, PAIRS))


def real_obstruction(v, samples=4001):
    """For a real AME state: (k, q) where k = dim of real first-order deformations modulo local orthogonal maps and
    q = min over unit deformation directions w of |P_coker (W W^T)_{cuts}|, the second-order obstruction.
    A real curve of AME states leaving the O(6)^4 orbit tangentially to w needs P_coker(W W^T) = 0, so q > 0
    means no real family through v (modulo local orthogonal maps).  k <= 2 is scanned on a fine circle."""
    v = np.asarray(v, float) / np.linalg.norm(v)
    J = np.vstack([_jac_block(v, A)[: D ** 4, 0::2].toarray() for A in [(0, 1), (0, 2), (0, 3)]]).real
    U, s, Vt = np.linalg.svd(J, full_matrices=True)
    r = int(np.sum(s > 1e-9 * s[0]))
    null, coker = Vt[r:].T, U[:, r:]
    t = v.reshape((D,) * N)
    gens = []
    for j in range(N):
        for a, b in itertools.combinations(range(D), 2):
            X = np.zeros((D, D)); X[a, b], X[b, a] = 1, -1
            gens.append(np.moveaxis(np.tensordot(X, t, axes=([1], [j])), 0, j).ravel())
    T = np.array(gens).T
    Ut, st, _ = np.linalg.svd(T, full_matrices=False)
    Qo = Ut[:, st > 1e-9 * st[0]]
    P = null - Qo @ (Qo.T @ null)
    Up, sp, _ = np.linalg.svd(P, full_matrices=False)
    W = Up[:, sp > 1e-6]
    k = W.shape[1]
    if k == 0:
        return 0, None
    from .core import CUTS, cut_matrix

    def Q(w):
        return coker.T @ np.concatenate([(cut_matrix(w, c) @ cut_matrix(w, c).T).ravel() for c in CUTS])
    if k == 1:
        return 1, float(np.linalg.norm(Q(W[:, 0])))
    if k == 2:
        phis = np.linspace(0, np.pi, samples)
        return 2, float(min(np.linalg.norm(Q(np.cos(p) * W[:, 0] + np.sin(p) * W[:, 1])) for p in phis))
    rng = np.random.default_rng(0)
    xs = rng.normal(size=(samples, k)); xs /= np.linalg.norm(xs, axis=1, keepdims=True)
    return k, float(min(np.linalg.norm(Q(W @ x)) for x in xs))
