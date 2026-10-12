"""Write paper/figures/amplitudes72.pdf: the three 2-unitaries of the uniform support-72 state Psi_0(0, pi/2).

Companion of figures/amplitudes.pdf (the real state psi, support 60).  At delta = 0, alpha = 0, beta = pi/2 the entries of
sqrt72 Psi are eps_x = +-1 on E, +1 on P and +i on P'.  Top row: the 36 x 36 reshapings for 12|34, 13|24, 14|23 in
natural order.  Bottom row: rows and columns permuted into the 18 blocks of size 2 x 2 -- first the six P/P' blocks,
with P on the main diagonal (they become the twelve 1 x 1 blocks of psi at delta = 1), then the twelve E blocks in
the order used for psi.  Colours: -1 and +1 as in amplitudes.pdf; +i amber (dataviz validate_palette.js,
--pairs all, light surface #fcfcfb: CVD dE 15.3, normal-vision dE 20.8; contrast 2.1:1, legend provides relief).

  python paper/make_amplitudes72.py [out.png]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))
from ame46real import core, family as fam  # noqa: E402

SURFACE, NEUTRAL, GRID = "#fcfcfb", "#f0efec", "#d9d8d4"
INK, INK2 = "#0b0b0b", "#52514e"
COLOR = {"-1": "#3987e5", "+1": "#e34948", "+i": "#eda100"}
LABEL = {"-1": r"$-1$", "+1": r"$+1$", "+i": r"$+i$"}

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 8,
                     "axes.edgecolor": GRID, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})

v = fam.psi72(1, 0.0, np.zeros(9))
assert core.ame_error(v) < 1e-12
w = np.sqrt(72) * v
KEY = {1: "+1", -1: "-1", 1j: "+i"}
state = {}
for q, x in enumerate(fam.POINTS):
    z = complex(np.round(w[fam.FLAT[q]], 12))
    state[x] = (KEY[z], {"L": "E", "O": "P", "N": "P'"}[fam.TYPE[q]])
assert len(state) == 72
pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]


def matrix(A, B):
    return {(6 * x[A[0]] + x[A[1]], 6 * x[B[0]] + x[B[1]]): s for x, s in state.items()}


def block_order(M):
    """Row and column orders: the P/P' blocks first (P on the main diagonal), then the E blocks."""
    rows = {}
    for (r, c), (_, t) in M.items():
        rows.setdefault(r, []).append((c, t))
    rord, cord, kinds, done = [], [], [], set()
    for kind in ("PP'", "E"):
        for r in sorted(rows):
            if r in done or (kind == "E") != (rows[r][0][1] == "E"):
                continue
            cs = sorted(c for c, _ in rows[r])
            partner = next(q for q in rows if q != r and sorted(c for c, _ in rows[q]) == cs)
            assert all(t == "E" for _, t in rows[r] + rows[partner]) == (kind == "E")
            if kind == "PP'":                      # put the P entry of row r first, so P lies on the diagonal
                cs = sorted(cs, key=lambda c: M[(r, c)][1] != "P")
            done |= {r, partner}
            rord += [r, partner]; cord += cs; kinds.append(kind)
    assert len(rord) == len(set(rord)) == 36 and len(cord) == len(set(cord)) == 36
    return rord, cord, kinds


def draw(ax, M, rpos, cpos, backing=None):
    ax.set_facecolor(SURFACE)
    if backing:
        for k in range(len(backing)):
            ax.add_patch(Rectangle((2 * k, 2 * k), 2, 2, facecolor=NEUTRAL, edgecolor="none", zorder=0))
    else:
        for k in range(0, 37, 6):
            ax.axhline(k, color=GRID, lw=0.4, zorder=0); ax.axvline(k, color=GRID, lw=0.4, zorder=0)
    inset = 0.09
    for (r, c), (val, _) in M.items():
        ax.add_patch(Rectangle((cpos[c] + inset, rpos[r] + inset), 1 - 2 * inset, 1 - 2 * inset,
                               facecolor=COLOR[val], edgecolor="none", zorder=2))
    ax.set_xlim(0, 36); ax.set_ylim(36, 0); ax.set_aspect("equal")
    ax.tick_params(length=0, pad=1.5)


fig, axes = plt.subplots(2, 3, figsize=(6.4, 4.75))
counts = []
for col, (A, B) in enumerate(pairings):
    M = matrix(A, B)
    ident = {k: k for k in range(36)}
    ax = axes[0, col]
    draw(ax, M, ident, ident)
    ax.set_xticks([3 + 6 * k for k in range(6)]); ax.set_xticklabels([str(k) for k in range(6)])
    ax.set_yticks([3 + 6 * k for k in range(6)]); ax.set_yticklabels([str(k) for k in range(6)])
    ax.set_title(rf"$\{{{A[0]+1},{A[1]+1}\}}\,|\,\{{{B[0]+1},{B[1]+1}\}}$", fontsize=9, pad=4)
    ax.set_xlabel(rf"column $(x_{B[0]+1}, x_{B[1]+1})$, blocks of $x_{B[0]+1}$", fontsize=7, color=INK2)
    ax.set_ylabel(rf"row $(x_{A[0]+1}, x_{A[1]+1})$, blocks of $x_{A[0]+1}$", fontsize=7, color=INK2, labelpad=1)
    rord, cord, kinds = block_order(M)
    counts.append((kinds.count("PP'"), kinds.count("E")))
    ax = axes[1, col]
    draw(ax, M, {r: i for i, r in enumerate(rord)}, {c: i for i, c in enumerate(cord)}, backing=kinds)
    ax.set_xticks([6, 12, 24, 36]); ax.set_yticks([6, 12, 24, 36])
    ax.set_xticklabels(["", "12", "", "36"]); ax.set_yticklabels(["", "12", "", "36"])
    ax.set_xlabel("permuted columns", fontsize=7, color=INK2)
    ax.text(14.5, 4.5, "six $2\\times2$ blocks\non $P\\cup P'$", fontsize=6.5, color=INK2, va="center")
    ax.text(1.0, 29.5, "twelve signed\n" + r"$2\times2$ blocks on $E$", fontsize=6.5, color=INK2)
    if col == 0:
        ax.set_ylabel("permuted rows", fontsize=7, color=INK2)
assert counts == [(6, 12)] * 3, counts

fig.text(0.012, 0.735, "natural order", rotation=90, va="center", fontsize=8, color=INK2)
fig.text(0.012, 0.29, "block order", rotation=90, va="center", fontsize=8, color=INK2)
order = ("-1", "+1", "+i")
handles = [Rectangle((0, 0), 1, 1, facecolor=COLOR[k], edgecolor="none") for k in order]
handles.append(Rectangle((0, 0), 1, 1, facecolor=SURFACE, edgecolor=GRID, lw=0.6))
fig.legend(handles, [LABEL[k] for k in order] + [r"$0$"], loc="lower center", ncol=4,
           frameon=False, fontsize=8, title=r"entries of $\sqrt{72}\,\Psi_0(0,\pi/2)$", title_fontsize=8,
           handlelength=1.1, handleheight=1.1, columnspacing=1.6, bbox_to_anchor=(0.5, -0.005))
fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.17, wspace=0.36, hspace=0.38)
fig.savefig(os.path.join(HERE, "figures", "amplitudes72.pdf"))
if len(sys.argv) > 1:
    fig.savefig(sys.argv[1], dpi=200)
print("blocks per pairing (P/P', E):", counts)
