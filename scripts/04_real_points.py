"""Real AME states on the support of psi: count, symmetry orbits, rigidity.

Exact:
  (1) real sign patterns on S60 = L u O solving the block congruences over GF(2): 2^32; modulo diagonal local
      phases they give exactly 2048 distinct states (the 2^11 choices theta_j in {0, pi} reach only 256 of them);
  (2) 512 of the 2048 lie on the gluing locus C (are limits of the support-72 family);
  (3) the support automorphism group Aut(S60) (party permutations x level relabellings) has order 192 and 16 orbits
      on the 2048 real states: 8 of size 64 and 8 of size 192; psi lies in an orbit of size 64.
Numerical:
  (4) real AME states on this support are isolated modulo local orthogonal maps: 14 of the 16 orbits (psi's
      included) have no real first-order deformation; at the other 2 the 2-dim space of real first-order
      deformations is obstructed at second order (|P_coker(W W^T)| >= 0.35 on every unit direction).
"""
from _common import check, done
from ame46real import realpoints as rp
from ame46real.jacobian import real_obstruction

signs, kdim, gdim = rp.real_points()
check("(1) 2^32 real sign patterns, gauge part 2^21 -> 2048 real states", kdim == 32 and gdim == 21 and len(signs) == 2048)
check("(1) theta_j in {0, pi} reach only 256 of them", len(rp.theta_points()) == 256)
check("(2) 512 real states on the gluing locus C", sum(rp.on_C(x) for x in signs) == 512)
perms = rp.automorphisms()
orbs = rp.orbits(signs, perms)
sizes = sorted(len(o) for o in orbs)
check("(3) |Aut(S60)| = 192", len(perms) == 192)
check("(3) 16 orbits: 8 x 64 + 8 x 192", sizes == [64] * 8 + [192] * 8)
check("(3) psi lies in an orbit of size 64", len(next(o for o in orbs if 0 in o)) == 64)
res = [real_obstruction(rp.state(signs[o[0]])) for o in orbs]
k0 = next(i for i, o in enumerate(orbs) if 0 in o)
print("    real first-order deformations per orbit:", [k for k, _ in res])
check("(4) psi: no real first-order deformation", res[k0][0] == 0)
check("(4) every orbit: no real deformation, or all of them obstructed at second order (numerical)",
      all(k == 0 or q > 1e-6 for k, q in res),
      "min obstruction " + ", ".join(f"{q:.3f}" for k, q in res if k) + " on orbits with " + str(sorted({k for k, _ in res if k})) + " directions")
done()
