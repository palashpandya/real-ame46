"""Real AME states on the two support-60 faces: counts, symmetry orbits, rigidity.

Face "+": support L u O (psi's face).  Face "-": support L u N (the delta = -1 end of the family; psi_minus).
Exact:
  (1) on each face, the real sign patterns solving the block congruences over GF(2) number 2^32; modulo diagonal local
      phases they give exactly 2048 distinct states (on face "+", the 2^11 choices theta_j in {0, pi} reach only 256);
  (2) 512 of them lie on the gluing locus (are limits of the support-72 family);
  (3) the support automorphisms (party permutations x level relabellings preserving support and officers) form a
      group of order 192 on each face, with 16 orbits on face "+" (8 x 64, 8 x 192) and 24 on face "-"
      (4 x 32, 20 x 96); the two faces are not related by any such monomial map.
Numerical:
  (4) the real states are isolated modulo local orthogonal maps: at every orbit representative there is either no
      real first-order deformation, or all of them are obstructed at second order.
"""
from _common import check, done
from ame46real import realpoints as rp
from ame46real.jacobian import real_obstruction

expected = {"+": (16, [64] * 8 + [192] * 8), "-": (24, [32] * 4 + [96] * 20)}
for s in "+-":
    F = rp.Face(s)
    signs, kdim, gdim = F.real_points()
    check(f"face {s}: (1) 2^32 real sign patterns, gauge part 2^21 -> 2048 real states",
          kdim == 32 and gdim == 21 and len(signs) == 2048)
    if s == "+":
        check("face +: (1) theta_j in {0, pi} reach only 256 of them", len(rp.theta_points()) == 256)
    check(f"face {s}: (2) 512 real states on the gluing locus", sum(F.on_C(x) for x in signs) == 512)
    perms = F.automorphisms()
    orbs = F.orbits(signs, perms)
    check(f"face {s}: (3) automorphism group of order 192 with {expected[s][0]} orbits",
          len(perms) == 192 and sorted(len(o) for o in orbs) == expected[s][1])
    res = [real_obstruction(F.state(signs[o[0]])) for o in orbs]
    print(f"    face {s}: real first-order deformations per orbit: {[k for k, _ in res]}")
    check(f"face {s}: (4) every orbit: no real deformation, or all obstructed at second order (numerical)",
          all(k == 0 or q > 1e-6 for k, q in res),
          "min obstructions " + (", ".join(f"{q:.3f}" for k, q in res if k) or "-"))
check("(3) the two faces are not monomially equivalent", not rp.monomially_equivalent(rp.Face("+"), rp.Face("-")))
done()
