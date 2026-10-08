"""The 2048 real states form exactly 8 classes under local unitaries and party permutations.

Lower bound (numerical; made rigorous in 07_certify_classes.py): the spectrum of 72 J^T J, J the Jacobian of
psi -> (rho_A) over all six 2-party marginals, is invariant under local unitaries (they act orthogonally on both
sides) and party permutations.  It takes exactly 8 distinct values on the 16 Aut(S60) orbits.
Upper bound (exact): the two orbits with equal spectrum are related by H_12^{(x)4} (the 2x2 Hadamard on levels
{1,2} of every party) followed by a monomial map.  Hence exactly 8 classes.  The representatives in
data/real_classes.json (psi with at most two officer signs flipped) are checked against the class table.
"""
import itertools

import numpy as np

from _common import check, done
from ame46real import core, realpoints as rp
from ame46real.jacobian import marginal_spectrum

signs, _, _ = rp.real_points()
orbs = rp.orbits(signs, rp.automorphisms())
orbit_of = {rp.key(signs[i]): k for k, o in enumerate(orbs) for i in o}
spec = [marginal_spectrum(rp.state(signs[o[0]])) for o in orbs]
dist = np.array([[np.max(np.abs(a - b)) for b in spec] for a in spec])
groups = []
for k in range(len(orbs)):
    for g in groups:
        if dist[k, g[0]] < 1e-9:
            g.append(k); break
    else:
        groups.append([k])
sep = min(dist[a, b] for g, h in itertools.combinations(groups, 2) for a in g for b in h)
noise = max(dist[a, b] for g in groups for a, b in itertools.combinations(g, 2))
check("lower bound: 8 distinct spectra on the 16 orbits, each shared by exactly 2 orbits",
      len(groups) == 8 and all(len(g) == 2 for g in groups), f"separation {sep:.2e}, noise {noise:.1e}")

H = np.eye(6); H[1:3, 1:3] = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
merged = 0
for g in groups:
    v = rp.state(signs[orbs[g[0]][0]]).reshape((6,) * 4)
    w = np.einsum("ai,bj,ck,dl,ijkl->abcd", H, H, H, H, v).ravel()
    x = rp.to_face_signs(w)
    merged += x is not None and orbit_of[rp.key(x)] == g[1]
check("upper bound: H_12^(x)4 + a monomial map joins the two orbits of each pair -> exactly 8 classes", merged == 8)

classes = core.load_json("real_classes.json")["classes"]
psi = core.psi()
print("\nclass  representative (officer signs flipped)  states  orbits      on C  N(41/8)")
ok = True
seen = set()
for c in classes:
    v = core.flip(psi, c["flip"])
    x = (v[rp.FLAT] < 0).astype(np.int64)
    k = orbit_of[rp.key(x)]
    g = next(i for i, gr in enumerate(groups) if k in gr)
    seen.add(g)
    members = [i for o in groups[g] for i in orbs[o]]
    nC = sum(rp.on_C(signs[i]) for i in members)
    cnt = int(np.sum(marginal_spectrum(v) < 41 / 8))
    ok &= (core.ame_error(v) < 1e-14 and len(members) == c["states"] and sorted(len(orbs[o]) for o in groups[g]) == c["aut_orbit_sizes"]
           and (nC == len(members)) == c["on_gluing_locus"] and (nC in (0, len(members))))
    print(f"{c['cls']:>5}  {' '.join(c['flip']) or '-':<38} {len(members):>6}  {'+'.join(str(len(orbs[o])) for o in groups[g]):<10} "
          f"{'yes' if nC else 'no':>4}  {cnt:>7}")
check("the 8 representatives are AME, lie in 8 different classes, and match the class table", ok and len(seen) == 8)
done()
