"""The -1 eigenspace of u^{(x)4}, u = diag(i,-1,-1,-i,-i,-i): random-start AME search and support minimisation over
local frames in the commutant U(1) x U(2) x U(3) of u (levels {0}, {1,2}, {3,4,5})."""
import itertools

import numpy as np
from scipy.optimize import minimize

from .core import CUTS, D, DIM, N, index

E = np.array([1, 2, 2, 3, 3, 3])
SEC = np.array([index(x) for x in itertools.product(range(D), repeat=N) if E[list(x)].sum() % 4 == 2])
MASK = np.zeros(DIM, bool); MASK[SEC] = True
PERMS = [c[0] + c[1] for c in CUTS]
INV = [np.argsort(p) for p in PERMS]


def embed(x):
    n = len(SEC); v = np.zeros(DIM, complex); v[SEC] = x[:n] + 1j * x[n:]; return v


def loss_grad(x):
    """f = sum over cuts of vdot(Delta, Delta), Delta = M M^dag - I/36 (normalised state); gradient in (Re, Im)."""
    v = embed(x); nv = np.linalg.norm(v); psi = v / nv
    t = psi.reshape((D,) * N); f = 0.0; G = np.zeros((D,) * N, complex)
    for p, ip in zip(PERMS, INV):
        M = np.transpose(t, p).reshape(36, 36)
        Dl = M @ M.conj().T - np.eye(36) / 36
        f += np.vdot(Dl, Dl).real
        G += np.transpose((4 * Dl @ M).reshape((D,) * N), ip)
    G = G.ravel(); G = (G - np.vdot(psi, G).real * psi) / nv
    return f, np.concatenate([G[SEC].real, G[SEC].imag])


def polar_steps(psi, steps=300):
    for _ in range(steps):
        for p, ip in zip(PERMS, INV):
            M = np.transpose(psi.reshape((D,) * N), p).reshape(36, 36)
            U, _, Vh = np.linalg.svd(M)
            psi = np.transpose((U @ Vh / 6).reshape((D,) * N), ip).ravel()
            psi[~MASK] = 0; psi /= np.linalg.norm(psi)
    return psi


def polish(x, stages=6):
    for _ in range(stages):
        f = loss_grad(x)[0]
        if f < 1e-28:
            break
        sc = min(max(1 / max(f, 1e-30), 1.0), 1e24)          # keep tiny losses well scaled for L-BFGS-B
        x = minimize(lambda y: tuple(np.multiply(sc, a) for a in loss_grad(y)), x, jac=True, method="L-BFGS-B",
                     options=dict(maxiter=5000, ftol=1e-15, gtol=1e-12)).x
    return x


def restart(rng):
    """One restart: Gaussian start in the sector, polar steps, rescaled L-BFGS.  Returns the normalised state."""
    n = len(SEC)
    psi = np.zeros(DIM, complex); psi[SEC] = rng.normal(size=n) + 1j * rng.normal(size=n)
    psi = polar_steps(psi / np.linalg.norm(psi))
    v = embed(polish(np.concatenate([psi[SEC].real, psi[SEC].imag])))
    return v / np.linalg.norm(v)


# ---------- support minimisation over commutant frames ------------------------------------------------------
BLOCKS = [[0], [1, 2], [3, 4, 5]]


def commutant_basis():
    out = []
    for b in BLOCKS:
        for i, a in enumerate(b):
            for c in b[i:]:
                if a == c:
                    H = np.zeros((D, D), complex); H[a, a] = 1; out.append(H)
                else:
                    H = np.zeros((D, D), complex); H[a, c] = H[c, a] = 1; out.append(H)
                    H = np.zeros((D, D), complex); H[a, c], H[c, a] = 1j, -1j; out.append(H)
    return np.array(out)                                      # 14 Hermitian generators


HB = commutant_basis()
K_ = len(HB)


def apply(Us, T):
    """(U_1 (x) U_2 (x) U_3 (x) U_4) T for a 6x6x6x6 tensor T (None = identity on that party)."""
    for j, U in enumerate(Us):
        if U is not None:
            T = np.moveaxis(np.tensordot(U, T, axes=([1], [j])), 0, j)
    return T


def _exp_and_dexp(h):
    """U = expm(i H) for H = sum_a h_a HB[a], and dU/dh_a for all a (closed form from the eigendecomposition)."""
    H = np.tensordot(h, HB, axes=1)
    lam, Q = np.linalg.eigh(H)
    el = np.exp(1j * lam)
    dl = lam[:, None] - lam[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        Gm = np.where(np.abs(dl) > 1e-12, (el[:, None] - el[None, :]) / np.where(np.abs(dl) > 1e-12, dl, 1),
                      1j * el[:, None] * np.ones_like(dl))
    U = (Q * el) @ Q.conj().T
    Bt = np.einsum("ik,akl,lj->aij", Q.conj().T, HB, Q)                  # Q^dag B_a Q
    dU = np.einsum("ik,akl,lj->aij", Q, Gm[None] * Bt, Q.conj().T)
    return U, dU


def l1_and_grad(p, T, eps):
    """f(p) = sum_x sqrt(|v_x|^2 + eps^2), v = (U_1 (x) ... (x) U_4) T, U_j = expm(i sum_a p_ja H_a)."""
    ed = [_exp_and_dexp(p[j * K_:(j + 1) * K_]) for j in range(N)]
    Us = [u for u, _ in ed]
    V = apply(Us, T)
    f = np.sum(np.sqrt(np.abs(V) ** 2 + eps ** 2))
    g = V / np.sqrt(np.abs(V) ** 2 + eps ** 2)                # df = Re sum conj(g) dV
    grad = np.zeros_like(p)
    for j in range(N):
        W = apply([Us[k] if k != j else None for k in range(N)], T)      # site j untouched
        Mj = np.tensordot(np.conj(np.moveaxis(g, j, 0)), np.moveaxis(W, j, 0), axes=([1, 2, 3], [1, 2, 3]))
        grad[j * K_:(j + 1) * K_] = np.real(np.einsum("aij,ij->a", ed[j][1], Mj))
    return f, grad


def frame(p):
    return [_exp_and_dexp(p[j * K_:(j + 1) * K_])[0] for j in range(N)]


def sparsify(v, rng, restarts=4, eps_schedule=(1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5), maxiter=2000):
    """Minimise the smoothed L1 norm over commutant frames.  Returns (support size at 1e-6, rotated state)."""
    T = np.asarray(v, complex).reshape((D,) * N)
    best = (10 ** 9, None)
    for r in range(restarts):
        p = rng.normal(size=N * K_) if r else np.zeros(N * K_)
        for eps in eps_schedule:
            p = minimize(l1_and_grad, p, args=(T, eps), jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter)).x
        w = apply(frame(p), T).ravel()
        supp = int(np.sum(np.abs(w) > 1e-6))
        if supp < best[0]:
            best = (supp, w)
    return best
