"""Reproduce the second discovery step: from a sector hit to the support-72 family (numerical).

1. Search the -1 eigenspace of u^{(x)4} (as in 08) until a hit with first-order deformation dimension 10 appears.
2. Minimise a smoothed L1 norm of the coefficients over local frames in the commutant U(1) x U(2) x U(3) of u (these
   frames keep the state in the eigenspace).
3. Keep the 72 largest coefficients and polish by Gauss-Newton on that fixed support to an exact AME state.
4. Check that this support is the array L u O u N up to a relabelling of parties and levels.  By Proposition 2(a),
   the state is then psi_sigma(delta, theta) up to diagonal local unitaries; delta is read off from the moduli.
Stochastic: the L1 minimisation reaches support 72 in a fraction of its restarts.

  python scripts/10_sparsify.py [seed=0]
"""
import itertools
import sys
import time

import numpy as np

from _common import check, done
from ame46real import core, family as fam, sector as sc
from ame46real.jacobian import _jac_block, deformation_dim
from ame46real.realpoints import relabelings

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
rng = np.random.default_rng(seed)
t0 = time.time()
for r in range(200):
    v = sc.restart(rng)
    if core.ame_error(v) < 1e-10 and deformation_dim(v)[0] == 10:
        break
print(f"sector hit with deformation dim 10 at restart {r}: support {len(core.support(v))} ({time.time() - t0:.0f} s)",
      flush=True)

found = None
for attempt in range(4):
    supp, w = sc.sparsify(v, rng, restarts=4)
    a = np.sort(np.abs(w))[::-1]
    print(f"L1 attempt {attempt}: support {supp} at 1e-6; |coefficients| 72nd / 73rd: {a[71]:.2e} / {a[72]:.2e} "
          f"({time.time() - t0:.0f} s)", flush=True)
    if supp == 72:
        found = w
        break
check("commutant frame with support 72 found", found is not None)
if found is None:
    done()

# Gauss-Newton on the fixed support
S = np.argsort(-np.abs(found))[:72]
x = found[S].copy()


def residual(x):
    u = np.zeros(core.DIM, complex); u[S] = x
    out = []
    for c in core.CUTS:
        M = core.cut_matrix(u, c); R = M @ M.conj().T - np.eye(36) / 36
        out.append(np.concatenate([R.real.ravel(), R.imag.ravel()]))
    return np.concatenate(out), u


cols = np.ravel(np.column_stack([2 * S, 2 * S + 1]))
for it in range(30):
    F, u = residual(x)
    if np.linalg.norm(F) < 1e-15:
        break
    J = np.vstack([_jac_block(u, c[0]).toarray() for c in core.CUTS])[:, cols]
    dx = np.linalg.lstsq(J, -F, rcond=None)[0]
    x = x + dx[0::2] + 1j * dx[1::2]
F, u = residual(x)
check("Gauss-Newton gives an exact AME state with exactly these 72 coefficients",
      core.ame_error(u) < 1e-14 and np.min(np.abs(x)) > 1e-3, f"AME error {core.ame_error(u):.1e}, "
      f"moved {np.linalg.norm(u - found):.1e} from the L1 frame")

Y = core.support(u)
target = set(fam.POINTS)
maps = [(sp, pi) for sp in itertools.permutations(range(4)) for pi in relabelings(Y, target, sp, first_only=True)]
check("the support is the array L u O u N up to a relabelling of parties and levels", len(maps) > 0)
if maps:
    sp, pi = maps[0]
    img = {tuple(pi[j][y[sp[j]]] for j in range(4)): u[core.index(y)] for y in Y}
    s = {t: sorted(72 * abs(img[fam.POINTS[q]]) ** 2 for q in fam.S72 if fam.TYPE[q] == t) for t in "LON"}
    delta = np.mean(s["O"]) - 1
    ok = (np.allclose(s["L"], 1, atol=1e-9) and np.allclose(s["O"], 1 + delta, atol=1e-9)
          and np.allclose(s["N"], 1 - delta, atol=1e-9))
    check("moduli: 1 on L, 1 + delta on O, 1 - delta on N (a member of the support-72 family)", ok,
          f"delta = {delta:+.4f}")
done()
