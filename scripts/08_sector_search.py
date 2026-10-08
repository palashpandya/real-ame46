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
import sys
import time

import numpy as np

import _common  # noqa: F401  (path setup)
from ame46real import core, sector as sc
from ame46real.jacobian import deformation_dim

R = int(sys.argv[1]) if len(sys.argv) > 1 else 30
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)


t0, hits = time.time(), []
for r in range(R):
    v = sc.restart(rng)
    err = core.ame_error(v)
    if err < 1e-10:
        dim = deformation_dim(v)[0]
        hits.append(dim)
        print(f"restart {r:3d}: HIT  AME error {err:.1e}, support {len(core.support(v))}, deformation dim {dim}", flush=True)
    else:
        print(f"restart {r:3d}: AME error {err:.1e}", flush=True)
print(f"\n{len(hits)} hits in {R} restarts ({time.time() - t0:.0f} s); deformation dims {sorted(hits)}")
print("dimension 10 = the support-72 family; other values are AME states of other families in the same sector")
