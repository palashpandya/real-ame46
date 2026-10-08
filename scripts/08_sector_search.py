"""Reproduce the discovery: search for AME(4,6) states in the -1 eigenspace of u^{(x)4}, u = diag(i,-1,-1,-i,-i,-i).

Only N = 4 and d = 6 enter.  The original search used the conjugated symmetry V u V^dag with Haar-random V; the loss
is LU-invariant, so searching in the diagonal frame is equivalent.  Each restart:
  Gaussian start in the 360-dim sector -> 300 sweeps of polar steps (replace each cut matrix by 1/6 times its unitary
  polar factor, project back to the sector) -> L-BFGS on sum_cuts |rho_A - I/36|^2 with rescaling (staged).
A restart is a hit when the AME error is < 1e-10.  For each hit the first-order deformation dimension is reported:
the original hit (support 348 in this frame) had dimension 10 and was brought to the support-72 form of the family
by optimising the local frame inside the commutant U(1) x U(2) x U(3) of u (see README).
Stochastic; typical: a few hits per 30 restarts, ~1 s per restart.

  python scripts/08_sector_search.py [restarts=30] [seed=0]
"""
import itertools
import sys
import time

import numpy as np
from scipy.optimize import minimize

import _common  # noqa: F401  (path setup)
from ame46real import core
from ame46real.jacobian import deformation_dim

R = int(sys.argv[1]) if len(sys.argv) > 1 else 30
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
e = np.array([1, 2, 2, 3, 3, 3])
sec = np.array([core.index(x) for x in itertools.product(range(6), repeat=4) if e[list(x)].sum() % 4 == 2])
n = len(sec)
mask = np.zeros(core.DIM, bool); mask[sec] = True
perms = [c[0] + c[1] for c in core.CUTS]
inv = [np.argsort(p) for p in perms]


def embed(x):
    v = np.zeros(core.DIM, complex); v[sec] = x[:n] + 1j * x[n:]; return v


def loss_grad(x):
    """f = sum_cuts vdot(Delta, Delta), Delta = M M^dag - I/36 for the normalised state; gradient in (Re, Im)."""
    v = embed(x); nv = np.linalg.norm(v); psi = v / nv
    t = psi.reshape((6,) * 4); f = 0.0; G = np.zeros((6,) * 4, complex)
    for p, ip in zip(perms, inv):
        M = np.transpose(t, p).reshape(36, 36)
        Dl = M @ M.conj().T - np.eye(36) / 36
        f += np.vdot(Dl, Dl).real
        G += np.transpose((4 * Dl @ M).reshape((6,) * 4), ip)
    G = G.ravel(); G = (G - np.vdot(psi, G).real * psi) / nv
    return f, np.concatenate([G[sec].real, G[sec].imag])


def polar_steps(psi, steps=300):
    for _ in range(steps):
        for p, ip in zip(perms, inv):
            M = np.transpose(psi.reshape((6,) * 4), p).reshape(36, 36)
            U, _, Vh = np.linalg.svd(M)
            psi = np.transpose((U @ Vh / 6).reshape((6,) * 4), ip).ravel()
            psi[~mask] = 0; psi /= np.linalg.norm(psi)
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


t0, hits = time.time(), []
for r in range(R):
    psi = np.zeros(core.DIM, complex); psi[sec] = rng.normal(size=n) + 1j * rng.normal(size=n)
    psi = polar_steps(psi / np.linalg.norm(psi))
    v = embed(polish(np.concatenate([psi[sec].real, psi[sec].imag])))
    v /= np.linalg.norm(v)
    err = core.ame_error(v)
    if err < 1e-10:
        dim = deformation_dim(v)[0]
        hits.append(dim)
        print(f"restart {r:3d}: HIT  AME error {err:.1e}, support {len(core.support(v))}, deformation dim {dim}", flush=True)
    else:
        print(f"restart {r:3d}: AME error {err:.1e}", flush=True)
print(f"\n{len(hits)} hits in {R} restarts ({time.time() - t0:.0f} s); deformation dims {sorted(hits)}")
print("dimension 10 = the support-72 family; other values are AME states of other families in the same sector")
