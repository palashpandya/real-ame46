"""Which symmetry sectors were searched: the necessary conditions that select them (inputs N = 4, d = 6 only).

A local symmetry U^{(x)4}, U_i = V_i diag(w^{e_i}) V_i^dag with w = exp(2 pi i / n), is described up to local frames by
the multisets of exponents e_i (one per party) and the eigenvalue w^s.  An AME state in the w^s eigenspace must pass:
  1. spectral test: for every set A of parties with |A| <= 2, the eigenvalue multiset of U_A must be contained in
     w^s times the conjugate eigenvalue multiset of U_{A^c}, with equality when |A| = 2 (rho_A is maximally mixed and
     U_A (x) U_{A^c} psi = w^s psi);
  2. population LP: in the eigenbasis, p_x = |psi_x|^2 is a probability vector on the allowed labels whose 2-party
     marginals are all uniform (the diagonal of rho_A = I/36), so this linear program must be feasible.
Classes are identified up to party permutations, per-party phases and units of Z_n; an order-n class whose exponents
share a factor g is the same symmetry as the order n/g class.  Expected (orders n <= 6): 204 distinct classes pass the
spectral test, 121 of them the LP, among them the order-4 class (1,2,2,3,3,3)^4 with s = 2, of dimension 360.

  python scripts/09_sector_filters.py [max_order=6]
"""
import itertools
import math
import sys
import time

import numpy as np
from scipy.optimize import linprog

import _common  # noqa: F401
from _common import check, done

N, d = 4, 6
LABELS = np.array(list(itertools.product(range(d), repeat=N)))
PAIRS = list(itertools.combinations(range(N), 2))


def conv(a, b, n):
    out = np.zeros(n, dtype=np.int64)
    for i in np.flatnonzero(a):
        out += a[i] * np.roll(b, i)
    return out


def hist(ms, sites, n):
    h = np.zeros(n, dtype=np.int64); h[0] = 1
    for i in sites:
        h = conv(h, ms[i], n)
    return h


def spectral_ok(ms, s, n):
    for a in range(N // 2, 0, -1):
        for A in itertools.combinations(range(N), a):
            Ac = [i for i in range(N) if i not in A]
            hA, r = hist(ms, A, n), hist(ms, Ac, n)[(s - np.arange(n)) % n]
            if (2 * a == N and not np.array_equal(hA, r)) or (2 * a < N and np.any(hA > r)):
                return False
    return True


def site_reps(n):
    reps, seen = [], set()
    for c in itertools.product(range(d + 1), repeat=n):
        if sum(c) == d:
            key = min(tuple(np.roll(c, j)) for j in range(n))
            if key not in seen:
                seen.add(key); reps.append(np.array(key))
    return reps


def canon_site(m, n):
    rots = [tuple(int(x) for x in np.roll(m, j)) for j in range(n)]
    best = min(rots)
    ks = [j for j in range(n) if rots[j] == best]
    return best, ks[0], (ks[1] - ks[0] if len(ks) > 1 else n)


def class_key(ms, s, n):
    best = None
    for a in range(1, n):
        if math.gcd(a, n) != 1:
            continue
        tot, g, sites = 0, n, []
        for m in ms:
            ma = np.zeros(n, dtype=np.int64)
            for t in range(n):
                ma[(a * t) % n] += m[t]
            c, k, p = canon_site(ma, n)
            sites.append(c); tot += k; g = math.gcd(g, p)
        cand = (tuple(sorted(sites)), (a * s + tot) % g, g)
        best = cand if best is None or cand < best else best
    return best


def exps_of(m):
    return [t for t in range(len(m)) for _ in range(int(m[t]))]


def reduce_class(n, exps, s):
    e = [np.array(x) - x[0] for x in exps]
    s2 = (s - sum(x[0] for x in exps)) % n
    g = n
    for x in e:
        for v in x:
            g = math.gcd(g, int(v))
    g = math.gcd(g, s2)
    n2 = n // g
    return n2, class_key([np.bincount((x // g) % n2, minlength=n2) for x in e], s2 // g, n2)


def sector(exps, n, s):
    E = sum(np.asarray(e)[LABELS[:, i]] for i, e in enumerate(exps))
    return np.flatnonzero(E % n == s)


def lp_ok(S):
    rows = []
    for A in PAIRS:
        idx = LABELS[S][:, A[0]] * d + LABELS[S][:, A[1]]
        M = np.zeros((d * d, len(S))); M[idx, np.arange(len(S))] = 1; rows.append(M)
    res = linprog(np.zeros(len(S)), A_eq=np.vstack(rows), b_eq=np.full(len(PAIRS) * d * d, 1 / d ** 2),
                  bounds=(0, None), method="highs")
    return res.status == 0


def enumerate_n(n):
    reps = site_reps(n)
    w = np.exp(2j * np.pi / n)
    F = np.array([[abs(np.sum(m * w ** (k * np.arange(n)))) for k in range(n)] for m in reps])
    Z = np.isclose(F, 0)
    out = {}
    # prefilter from the balanced cuts: at each frequency the four character moduli agree, or >= 3 vanish
    for idx in itertools.combinations_with_replacement(range(len(reps)), N):
        f, z = F[list(idx)], Z[list(idx)]
        if not (np.all(np.isclose(f, f[0]), axis=0) | (z.sum(axis=0) >= 3)).all():
            continue
        ms = [reps[i] for i in idx]
        for s in range(n):
            if spectral_ok(ms, s, n):
                exps = [exps_of(m) for m in ms]
                out.setdefault(reduce_class(n, exps, s), (n, exps, s))
    return out


max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
allc, t0 = {}, time.time()
for n in range(1, max_n + 1):
    new = {k: v for k, v in enumerate_n(n).items() if k not in allc}
    for key, (n_, exps, s) in new.items():
        S = sector(exps, n_, s)
        allc[key] = dict(n=n_, exps=exps, s=s, dim=len(S), lp=lp_ok(S))
    print(f"order {n}: {len(new)} new classes pass the spectral test ({time.time() - t0:.0f} s)", flush=True)
lp = [c for c in allc.values() if c["lp"]]
target = [c for c in lp if c["exps"] == [[1, 2, 2, 3, 3, 3]] * 4 and c["s"] == 2]
if max_n == 6:
    check("204 classes pass the spectral test, 121 of them the population LP", len(allc) == 204 and len(lp) == 121,
          f"{len(allc)}, {len(lp)}")
check("the order-4 class (1,2,2,3,3,3)^4, eigenvalue w^2 = -1, passes both, with sector dimension 360",
      len(target) == 1 and target[0]["dim"] == 360)
done()
