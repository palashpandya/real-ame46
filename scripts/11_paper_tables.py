"""The tables printed in the paper (paper/inputs/*.tex and Table I in paper/main.tex) agree with the data.

Checks: the displayed state, the sets O and N, the table of K, the branch table, the three block tables of
Appendix A (every 1x1 entry and every signed 2x2 block against psi), and the class table of Theorem 1.
"""
import os
import re

from _common import check, done
from ame46real import core, family as fam

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "paper")
read = lambda f: open(os.path.join(P, f)).read()

terms = {}
for row in re.findall(r"^\$a=(\d)\$ & (.*?)\\\\$", read("inputs/state_display.tex"), re.M):
    a, cells = row
    for b, cell in enumerate(cells.split(" & ")):
        for sign, root, cd in re.findall(r"([+-]?)(\\sqrt2\\,)?\\ket\{(\d\d)\}", cell):
            terms[a + str(b) + cd] = (sign or "+") + ("sqrt2" if root else "1")
check("displayed state (officer table) = data/psi_terms.json", terms == dict(core.load_json("psi_terms.json")["terms"]), f"{len(terms)} terms")

st = read("inputs/family_sets.tex")
O = sorted(core.lab(fam.POINTS[q]) for q in fam.S72 if fam.TYPE[q] == "O")
Nn = sorted(core.lab(fam.POINTS[q]) for q in fam.S72 if fam.TYPE[q] == "N")
check("sets O and N", sorted(re.findall(r"\d{4}", st.split("N&")[0])) == O and sorted(re.findall(r"\d{4}", st.split("N&")[1])) == Nn)

Kd = {}
for j, plus, minus in re.findall(r"\$\\theta_\{(\d)\}\$ & (.*?) & (.*?)\\\\", read("inputs/family_K.tex")):
    Kd.update({(x, int(j)): 1 for x in re.findall(r"(\d{4})", plus)})
    Kd.update({(x, int(j)): -1 for x in re.findall(r"(\d{4})", minus)})
check("table of K", Kd == {(core.lab(fam.POINTS[q]), j + 1): int(fam.K[q, j]) for q in fam.S72 for j in range(9) if fam.K[q, j]})

bt = read("inputs/family_branches.tex")
th = [("0",) * 9] + [tuple(m.split(",")) for m in re.findall(r"\\pi\(([\d,]+)\)", bt)]
check("branch table", {(s, t) for t in th for s in (1, -1)} == {(b["sigma"], tuple(b["theta_over_pi"])) for b in fam.F["branches"]})

psi_t = core.psi_terms()
for part, cut in zip(read("inputs/blocks.tex").split("\\paragraph")[1:], core.CUTS):
    ones = re.findall(r"\$(\d\d)\\to(\d\d)\$", part)
    ok = len(ones) == 12 and all(psi_t.get(core.point(tuple(map(int, r)), tuple(map(int, c)), cut)) == (0, 1) for r, c in ones)
    twos = re.findall(r"\$(\d\d),(\d\d)\$; \$(\d\d),(\d\d)\$; \$\\left\(\\begin\{smallmatrix\}([+-])&([+-])\\\\([+-])&([+-])\\end", part)
    ok &= len(twos) == 12
    for r1, r2, c1, c2, *sg in twos:
        vals = [psi_t.get(core.point(tuple(map(int, r)), tuple(map(int, c)), cut)) for r in (r1, r2) for c in (c1, c2)]
        ok &= all(v is not None and v[1] == 0 and v[0] == (1 if s == "+" else -1) for v, s in zip(vals, sg))
        ok &= (1 if sg[0] == "+" else -1) * (1 if sg[2] == "+" else -1) + (1 if sg[1] == "+" else -1) * (1 if sg[3] == "+" else -1) == 0
    check(f"Appendix A block table for {''.join(str(a + 1) for a in cut[0])}|{''.join(str(a + 1) for a in cut[1])}", ok)

rows = re.findall(r"^(\d+) & \$([+-])\$ & (.*?) & \$(\d+)\$ & (yes|no) & \$(\d+)\$ & \$(\d+)\$ & \$(\d+)\$\\\\", read("main.tex"), re.M)
data = {(c["cls"], c["face"], tuple(c["flip"]), c["states"], c["on_gluing_locus"], c["deformation_dim"], c["N_149_8"], c["N_3_2"])
        for c in core.load_json("real_classes.json")["classes"]}
tex = {(int(a), f, tuple(re.findall(r"\d{4}", fl)), int(s), y == "yes", int(dm), int(n1), int(n2)) for a, f, fl, s, y, dm, n1, n2 in rows}
check("Table I (classes) = data/real_classes.json", tex == data, f"{len(tex)} rows")
done()
