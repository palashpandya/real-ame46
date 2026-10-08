"""The real AME(4,6) state psi: exact verification in Q(sqrt2), block structure, officers.

  sqrt(72) psi has 48 coefficients +-1 and 12 coefficients sqrt2 (data/psi_terms.json).
Checks
  (1) <psi|psi> = 1 and rho_A = I/6^|A| EXACTLY for all 4 one-party and 6 two-party subsets A (Q(sqrt2) arithmetic).
  (2) In each of the three 2|2 reshapings, O = 6 M is (up to permutations) a direct sum of 12 blocks (1) and 12 blocks
      (1/sqrt2)[[s1, s2], [s3, s4]] with s_i = +-1 and s1 s3 + s2 s4 = 0; hence O is a real orthogonal 2-unitary.
  (3) Officers: in every reshaping, the 36 rows (and the 36 columns) are 12 basis vectors ("classical officers") and
      24 two-term states (|kl> +- |k'l'>)/sqrt2 with k != k', l != l' (Bell-type, Schmidt rank 2).
  (4) The 12 coefficients 1/6 saturate the bound |<abcd|psi>| <= 1/6 valid for every AME(4,6) state
      (|<abcd|psi>|^2 <= <ab|rho_12|ab> = 1/36): psi has the maximal possible overlap with product states.
"""
import itertools
from fractions import Fraction as Fr

from _common import check, done
from ame46real.core import CUTS, D, N, blocks, point, psi_terms

state = psi_terms()                                     # x -> (a, b): sqrt72 psi_x = a + b sqrt2
mul = lambda x, y: (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])
add = lambda x, y: (x[0] + y[0], x[1] + y[1])
check("support 60 = 48 x (+-1) + 12 x sqrt2",
      len(state) == 60 and sum(abs(a) == 1 and b == 0 for a, b in state.values()) == 48
      and sum(a == 0 and b == 1 for a, b in state.values()) == 12)

# (1) exact marginals (all entries scaled by 72)
norm = (0, 0)
for c in state.values():
    norm = add(norm, mul(c, c))
check("(1) <psi|psi> = 1 exactly", norm == (72, 0))
for r in (1, 2):
    for A in itertools.combinations(range(N), r):
        Ac = [i for i in range(N) if i not in A]
        rho = {}
        for s, cs in state.items():
            for t, ct in state.items():
                if all(s[i] == t[i] for i in Ac):
                    k = (tuple(s[i] for i in A), tuple(t[i] for i in A))
                    rho[k] = add(rho.get(k, (0, 0)), mul(cs, ct))
        target = Fr(72, D ** r)
        ok = all(v == ((target, 0) if k[0] == k[1] else (0, 0)) for k, v in rho.items())
        ok &= all((x, x) in rho for x in itertools.product(range(D), repeat=r))
        check(f"(1) rho_{''.join(str(a + 1) for a in A)} = I/{D ** r} exactly", ok)

# (2) block structure; the complete block lists are printed (they replace the table of the paper's first draft)
for cut in CUTS:
    nb1 = nb2 = 0
    name = "".join(str(a + 1) for a in cut[0]) + "|" + "".join(str(a + 1) for a in cut[1])
    ones, twos = [], []
    for R, C in blocks(state.keys(), cut):
        pt = lambda r, c: point(r, c, cut)
        lab2 = lambda t: "%d%d" % t
        if len(R) == len(C) == 1:
            nb1 += state[pt(R[0], C[0])] == (0, 1)
            ones.append(f"{lab2(R[0])}->{lab2(C[0])}")
        elif len(R) == len(C) == 2:
            M = [[state.get(pt(r, c)) for c in C] for r in R]
            if all(m is not None and m[1] == 0 and abs(m[0]) == 1 for row in M for m in row):
                sg = [[m[0] for m in row] for row in M]
                nb2 += sg[0][0] * sg[1][0] + sg[0][1] * sg[1][1] == 0
                sym = lambda v: "+" if v > 0 else "-"
                twos.append(f"rows {lab2(R[0])},{lab2(R[1])} cols {lab2(C[0])},{lab2(C[1])} "
                            f"[[{sym(sg[0][0])},{sym(sg[0][1])}],[{sym(sg[1][0])},{sym(sg[1][1])}]]")
    print(f"    cut {name}: 1x1 blocks (row->column, entry 1 of O = 6M): {', '.join(ones)}")
    print(f"    cut {name}: 2x2 blocks (1/sqrt2 x sign matrix):")
    for t in twos:
        print("        " + t)
    check(f"(2) cut {name}: 12 blocks (1) + 12 orthogonal signed 2x2 Hadamard blocks", nb1 == 12 and nb2 == 12,
          f"{nb1} + {nb2}")

# (3) officers
for cut in CUTS:
    A_, B_ = cut
    for side, (P_, Q_) in (("rows", (A_, B_)), ("columns", (B_, A_))):
        lines = {}
        for x in state:
            lines.setdefault((x[P_[0]], x[P_[1]]), []).append((x[Q_[0]], x[Q_[1]]))
        classical = sum(len(v) == 1 for v in lines.values())
        bell = sum(len(v) == 2 and v[0][0] != v[1][0] and v[0][1] != v[1][1] for v in lines.values())
        name = "".join(str(a + 1) for a in A_) + "|" + "".join(str(a + 1) for a in B_)
        check(f"(3) cut {name}, {side}: 12 classical officers + 24 Bell-type officers",
              len(lines) == 36 and classical == 12 and bell == 24, f"{classical} + {bell}")

# (4) product overlap
officers = sorted("".join(map(str, x)) for x, c in state.items() if c == (0, 1))
check("(4) 12 product states with overlap 1/6 (the AME maximum)", len(officers) == 12, " ".join(officers))
done()
