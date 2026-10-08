"""The 4096 real states of the two faces form exactly 20 classes under local unitaries and party permutations.

Upper bound (exact): states in one automorphism orbit are related by monomial maps; orbits of the same face are
further joined by H_12^{(x)4} (the 2x2 Hadamard on levels {1,2} of every party) followed by a monomial map.  This
leaves 8 classes on face "+" and 12 on face "-".
Lower bound (numerical; made rigorous in 07_certify_classes.py): the spectrum of 72 J^T J, J the Jacobian of
psi -> (rho_A) over all six 2-party marginals, is invariant under local unitaries and party permutations.  It is
constant on each of the 20 classes and takes 20 distinct values.
The representatives in data/real_classes.json are checked against the class table.
"""
import itertools

import numpy as np

from _common import check, done
from ame46real import core, realpoints as rp
from ame46real.jacobian import marginal_spectrum

H = np.eye(6); H[1:3, 1:3] = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
table = core.load_json("real_classes.json")["classes"]
class_spec = []
for s in "+-":
    F = rp.Face(s)
    signs, _, _ = F.real_points()
    orbs = F.orbits(signs, F.automorphisms())
    orbit_of = {F.key(signs[i]): k for k, o in enumerate(orbs) for i in o}
    parent = list(range(len(orbs)))

    def find(a):
        while parent[a] != a:
            a = parent[a]
        return a
    merges = 0
    for k, o in enumerate(orbs):
        w = np.einsum("ai,bj,ck,dl,ijkl->abcd", H, H, H, H, F.state(signs[o[0]]).reshape((6,) * 4)).ravel()
        x = F.to_face_signs(w)
        if x is not None:
            a, b = find(k), find(orbit_of[F.key(x)])
            merges += a != b
            if a != b:
                parent[a] = b
    classes = {}
    for k in range(len(orbs)):
        classes.setdefault(find(k), []).append(k)
    spec = [marginal_spectrum(F.state(signs[o[0]])) for o in orbs]
    const = max((np.max(np.abs(spec[a] - spec[b])) for c in classes.values() for a, b in itertools.combinations(c, 2)),
                default=0.0)
    n_expected = 8 if s == "+" else 12
    check(f"face {s}: upper bound: H_12^(x)4 merges {len(orbs)} orbits into {n_expected} classes",
          len(classes) == n_expected, f"{merges} merges")
    check(f"face {s}: the spectrum is constant on each class (consistency)", const < 1e-9, f"max deviation {const:.1e}")
    rows = [r for r in table if r["face"] == s]
    ok, seen = True, set()
    for r in rows:
        x = F.flips(r["flip"])
        c = find(orbit_of[F.key(x)])
        seen.add(c)
        members = [i for k in classes[c] for i in orbs[k]]
        sp = marginal_spectrum(F.state(x))
        nC = sum(F.on_C(signs[i]) for i in members)
        ok &= (core.ame_error(F.state(x)) < 1e-14 and len(members) == r["states"]
               and sorted(len(orbs[k]) for k in classes[c]) == r["aut_orbit_sizes"]
               and nC == (len(members) if r["on_gluing_locus"] else 0)
               and int(np.sum(sp < 149 / 8)) == r["N_149_8"] and int(np.sum(sp < 1.5)) == r["N_3_2"])
        class_spec.append((r["cls"], sp))
    check(f"face {s}: the {len(rows)} representatives are AME, in different classes, and match the class table",
          ok and len(seen) == len(rows) == n_expected)

spec_of = dict(class_spec)
pt = [[np.sum(spec_of[c] ** k) for c in range(3, 9)] for k in range(1, 9)]
check("power traces Tr (72 J^T J)^k, k <= 8, do not separate classes 3-8 (numerical)",
      all(np.ptp(t) <= 1e-9 * max(abs(t[0]), 1) for t in pt))
sep = min(np.max(np.abs(a - b)) for (_, a), (_, b) in itertools.combinations(class_spec, 2))
check("lower bound: the 20 class representatives have 20 distinct spectra (both faces together)", sep > 1e-6,
      f"min separation {sep:.2e}")
print("\nclass face  flipped        states  orbits     on C  dim  N(149/8) N(3/2)")
for r in table:
    print(f"{r['cls']:>5}  {r['face']}    {' '.join(r['flip']) or '-':<14} {r['states']:>5}   "
          f"{'+'.join(map(str, r['aut_orbit_sizes'])):<10} {'yes' if r['on_gluing_locus'] else 'no':>4} "
          f"{r['deformation_dim']:>4} {r['N_149_8']:>8} {r['N_3_2']:>6}")
done()
