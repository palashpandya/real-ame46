"""Entangling power: every real AME state gives a real orthogonal two-quhex gate of maximal entangling power.

For a gate U on C^d (x) C^d the entangling power [Zanardi, Zalka, Faoro, PRA 62, 030301 (2000)] is the mean linear
entropy it produces from Haar-random product inputs,
    e_p(U) = < 1 - tr rho_A^2 >,   rho_A = tr_B U|a><a| (x) |b><b| U^dag,
and [Zanardi, PRA 63, 040304 (2001)]
    e_p(U) = (d/(d+1))^2 [E(U) + E(US) - E(S)],   E(U) = 1 - tr (U^R U^R^dag)^2 / d^4,   E(US) = 1 - tr (U^G U^G^dag)^2 / d^4,
with S the swap and E(S) = 1 - 1/d^2.  Hence e_p(U) <= (d-1)/(d+1) = 5/7, with equality iff U is 2-unitary
[Rather, Aravinda, Lakshminarayan, PRL 125, 070501 (2020)].  We take U = 6 psi reshaped along 12|34 (rows x_1 x_2,
columns x_3 x_4) and report the normalised entangling power e_p / (5/7).  Checks (numerical):
  (1) psi, U along each of the three cuts: e_p = 5/7, both from Zanardi's formula and from the exact Haar average over
      product inputs (Weingarten: E[a_p a_r a*_p' a*_r'] = (d_pp' d_rr' + d_pr' d_rp') / (d(d+1)); no use of the
      formula), and from a Monte Carlo average.  Since O, O^R and O^G are orthogonal exactly (01), e_p = 5/7 exactly;
  (2) all 4096 real states of both faces, U along every one of the 24 orderings of the parties: e_p / (5/7) = 1;
  (3) the same for random members of both families and for the golden 2-unitary of Rather et al.;
  (4) controls: the two computations agree, and give e_p < 5/7, for a Haar-random U(36) and O(36); SWAP gives 0;
  (5) e_p is an average: on the 36 basis inputs |x_3 x_4> the gate 6 psi gives 12 product outputs (the classical
      officers, E = 0) and 24 two-term Bell-type outputs (E = 1/2).
"""
import itertools
import os

import numpy as np

from _common import check, done
from ame46real import core, family as fam
from ame46real import realpoints as rp

d = core.D
EP_MAX = (d - 1) / (d + 1)                                          # 5/7
E_S = 1 - 1 / d ** 2                                                # 35/36


def gate(v, order=(0, 1, 2, 3)):
    """U = d psi with rows (x_i x_j) and columns (x_k x_l), (i, j, k, l) = order."""
    v = np.asarray(v) / np.linalg.norm(v)
    return d * np.transpose(v.reshape((d,) * 4), order).reshape(d * d, d * d)


def realign(U):                                                     # <n m|U^R|a b> = <n a|U|m b>
    return U.reshape(d, d, d, d).transpose(0, 2, 1, 3).reshape(d * d, d * d)


def ptrans(U):                                                      # <m a|U^G|n b> = <n a|U|m b>
    return U.reshape(d, d, d, d).transpose(2, 1, 0, 3).reshape(d * d, d * d)


def op_ent(X):
    G = X @ X.conj().T
    return 1 - np.trace(G @ G).real / d ** 4


def ep_formula(U):
    """Zanardi's formula; returns (e_p, E(U), E(US))."""
    EU, EUS = op_ent(realign(U)), op_ent(ptrans(U))
    return (d / (d + 1)) ** 2 * (EU + EUS - E_S), EU, EUS


def ep_haar(U):
    """1 - E_{a,b} tr rho_A^2 from the second moments of Haar-random |a>, |b> (exact; independent of the formula)."""
    T = U.reshape(d, d, d, d)                                       # T[i, j, p, q] = <i j|U|p q>
    tot = 0.0
    for pa, ra in (("p", "r"), ("r", "p")):                         # (p', r') paired with (p, r) or (r, p)
        for qb, sb in (("q", "s"), ("s", "q")):                     # (q', s') paired with (q, s) or (s, q)
            tot += np.einsum(f"ijpq,kj{pa}{qb},klrs,il{ra}{sb}->", T, T.conj(), T, T.conj(), optimize=True).real
    return 1 - tot / (d * (d + 1)) ** 2


def lin_entropies(U, inputs):
    """Linear entropy of U|a>|b> for each row (a (x) b) of inputs."""
    out = (inputs @ U.T).reshape(-1, d, d)
    rho = out @ out.conj().transpose(0, 2, 1)
    return 1 - np.einsum("nij,nji->n", rho, rho).real


def ep_monte_carlo(U, n, rng):
    a = rng.normal(size=(n, d)) + 1j * rng.normal(size=(n, d))
    b = rng.normal(size=(n, d)) + 1j * rng.normal(size=(n, d))
    a /= np.linalg.norm(a, axis=1, keepdims=True)
    b /= np.linalg.norm(b, axis=1, keepdims=True)
    e = lin_entropies(U, np.einsum("np,nq->npq", a, b).reshape(n, d * d))
    return e.mean(), e.std() / np.sqrt(n)


def haar(n, rng, real=False):
    Z = rng.normal(size=(n, n)) if real else rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


print(f"d = {d}: E(S) = 35/36 = {E_S:.6f}, maximal entangling power (d-1)/(d+1) = 5/7 = {EP_MAX:.6f}\n")
rng = np.random.default_rng(13)
psi = core.psi()

# (1) psi along the three cuts
dev = 0.0
for A, B in core.CUTS:
    U = gate(psi, A + B)
    ep, EU, EUS = ep_formula(U)
    eh = ep_haar(U)
    dev = max(dev, abs(ep - EP_MAX), abs(eh - EP_MAX))
    print(f"    psi, cut {A[0] + 1}{A[1] + 1}|{B[0] + 1}{B[1] + 1}: E(U) = {EU:.12f}, E(US) = {EUS:.12f}, "
          f"e_p = {ep:.12f} (formula), {eh:.12f} (Haar average)")
check("(1) psi: e_p = 5/7 along all three cuts, from the formula and from the exact Haar average", dev < 1e-12,
      f"max deviation {dev:.1e}")
m, se = ep_monte_carlo(gate(psi), 100_000, rng)
check("(1) psi: Monte Carlo average over 10^5 random product inputs agrees with 5/7 (within 4 s.e.)",
      abs(m - EP_MAX) < 4 * se, f"{m:.5f} +- {se:.5f}, normalised {m / EP_MAX:.5f}")

# (2) all real states, all orderings of the parties
worst, count = 0.0, 0
for s in "+-":
    F = rp.Face(s)
    signs, _, _ = F.real_points()
    for x in signs:
        v = F.state(x)
        for order in itertools.permutations(range(4)):
            worst = max(worst, abs(ep_formula(gate(v, order))[0] / EP_MAX - 1))
            count += 1
check(f"(2) all 4096 real states x 24 orderings of the parties ({count} gates): normalised e_p = 1", worst < 1e-12,
      f"max |e_p / (5/7) - 1| = {worst:.1e}")
table = core.load_json("real_classes.json")["classes"]
faces = {s: rp.Face(s) for s in "+-"}
dev = max(abs(ep_haar(gate(faces[r["face"]].state(faces[r["face"]].flips(r["flip"])))) - EP_MAX) for r in table)
check("(2) the 20 class representatives: exact Haar average = 5/7", dev < 1e-12, f"max deviation {dev:.1e}")

# (3) families and the golden 2-unitary
others = {"psi_sigma(delta, theta), random": gate(fam.psi72(1, rng.uniform(-1, 1), rng.uniform(0, 7, 9))),
          "psi_face(theta), random": gate(fam.psi_face(rng.uniform(0, 7, 11))),
          "golden 2-unitary (Rather et al. 2022)": np.load(os.path.join(core.DATA, "reference", "golden_rbb22_U36.npy"))}
for name, U in others.items():
    ep, eh = ep_formula(U)[0], ep_haar(U)
    check(f"(3) {name}: normalised e_p = 1", abs(ep / EP_MAX - 1) < 1e-12 and abs(eh / EP_MAX - 1) < 1e-12,
          f"{ep / EP_MAX:.12f}")

# (4) controls
swap = np.eye(d * d)[[d * (i % d) + i // d for i in range(d * d)]]
controls = {"Haar-random U(36)": haar(36, rng), "Haar-random O(36)": haar(36, rng, real=True), "SWAP": swap,
            "identity": np.eye(d * d)}
for name, U in controls.items():
    ep, EU, EUS = ep_formula(U)
    eh = ep_haar(U)
    print(f"    {name}: E(U) = {EU:.6f}, E(US) = {EUS:.6f}, normalised e_p = {ep / EP_MAX:.6f} "
          f"(Haar average {eh / EP_MAX:.6f})")
    ok = abs(ep - eh) < 1e-12 and (ep < EP_MAX - 1e-3 if name.startswith("Haar") else abs(ep) < 1e-12)
    check(f"(4) {name}: formula = Haar average, " + ("e_p < 5/7" if name.startswith("Haar") else "e_p = 0"), ok)

# (5) basis inputs
e = lin_entropies(gate(psi), np.eye(d * d))
n0, nh = int(np.sum(np.abs(e) < 1e-12)), int(np.sum(np.abs(e - 0.5) < 1e-12))
check("(5) psi on the 36 basis inputs: 12 product outputs (E = 0) and 24 outputs with E = 1/2 (mean 1/3 < 5/7)",
      n0 == 12 and nh == 24, f"{n0} + {nh}")
done()
