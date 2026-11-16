"""Exact real isolation and effective local-stabiliser Lie algebra for the two faces.

All real sign patterns (2^32 per face) are covered by local signs, support
permutations, and 512 explicitly checked candidates. Unlike the diagonal-LU
classification, this reduction uses only real orthogonal maps. At each candidate
integer identities over Z[sqrt(2)] and finite-field rank lower bounds certify
isolation in an unrestricted real slice. See the proposition in paper/main.tex.

Run from the repository root, also with python -I. No sympy or numerical kernel
finder is used. --workers defaults to 4, with at most 16 total processes.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import itertools
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time

# -I removes the script directory; resolve only this repository's helper path.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import check, done

# Bound CPU use independently of the legacy runner's BLAS settings.
for _thread_variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_thread_variable] = '1'
import numpy as np
from ame46real import family as fam, realpoints as rp
from ame46real.lattice import gf2_nullspace, gf2_rank

DIM = 1296
CUTS = ((0, 1), (0, 2), (0, 3))
PRIMES = (2147483647, 2147483543)
IDX = np.arange(DIM).reshape((6,) * 4)
TRI = np.triu_indices(36)
TARGET = np.tile((2 * np.eye(36, dtype=np.int64))[TRI], 3)


def rank_mod_p(matrix, p):
    """Gaussian elimination, reducing every product before its next use.

    For p < 2^31, products and subtractions stay in signed int64: each
    product is <= (p-1)^2 < 2^62. There are no unreduced dot products.
    """
    assert 2 < p < 2**31
    a = np.array(matrix, dtype=np.int64) % p
    m, n = a.shape
    r = 0
    for c in range(n):
        if r == m:
            break
        piv = np.flatnonzero(a[r:, c])
        if not piv.size:
            continue
        i = r + piv[0]
        a[[r, i]] = a[[i, r]]
        a[r, c:] = a[r, c:] * pow(int(a[r, c]), p - 2, p) % p
        rows = r + 1 + np.flatnonzero(a[r + 1:, c])
        if rows.size:
            ix = np.ix_(rows, np.arange(c, n))
            a[ix] = (a[ix] - a[rows, c, None] * a[r, c:]) % p
        r += 1
    return r


def jac(v):
    """Integer derivative of the upper triangles of the three V V^T."""
    out = []
    for cut in CUTS:
        axes = cut + tuple(i for i in range(4) if i not in cut)
        ids = IDX.transpose(axes).reshape(36, 36)
        vcut = v[ids]
        block = np.zeros((666, DIM), dtype=np.int64)
        for row, (i, j) in enumerate(zip(*TRI)):
            block[row, ids[i]] += vcut[j]
            block[row, ids[j]] += vcut[i]
        out.append(block)
    return np.vstack(out)


def quad(v):
    out = []
    for cut in CUTS:
        axes = cut + tuple(i for i in range(4) if i not in cut)
        vcut = v[IDX.transpose(axes).reshape(36, 36)]
        out.append((vcut @ vcut.T)[TRI])
    return np.concatenate(out)


def tangent(v, symmetric=False):
    cols = []
    for j in range(4):
        pairs = (itertools.combinations_with_replacement(range(6), 2) if symmetric
                 else itertools.combinations(range(6), 2))
        for a, b in pairs:
            h = np.zeros((6, 6), dtype=np.int64)
            h[a, b] = 1
            h[b, a] = 1 if symmetric else -1
            cols.append(np.moveaxis(np.tensordot(h, v.reshape((6,) * 4),
                                                axes=([1], [j])), 0, j).ravel())
    return np.array(cols).T


def rational_kernel(matrix):
    """Small exact rational RREF, used only on 12 candidate support columns."""
    a = [[Fraction(int(x)) for x in row] for row in matrix if np.any(row)]
    r, piv = 0, []
    n = matrix.shape[1]
    for c in range(n):
        k = next((i for i in range(r, len(a)) if a[i][c]), None)
        if k is None:
            continue
        a[r], a[k] = a[k], a[r]
        z = a[r][c]
        a[r] = [x / z for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [x - z * y for x, y in zip(a[i], a[r])]
        piv.append(c)
        r += 1
    out = []
    for c in sorted(set(range(n)) - set(piv)):
        w = [Fraction(0)] * n
        w[c] = Fraction(1)
        for i, p in enumerate(piv):
            w[p] = -a[i][c]
        # Clear denominators without sympy, using Python's exact integers.
        import math
        den = math.lcm(*(x.denominator for x in w))
        out.append(np.array([int(x * den) for x in w], dtype=np.int64))
    return out


def require(ok, message):
    if not ok:
        raise ValueError(message)


def certify(task):
    tag, oi, coset, x, flat, officer, other_cols, roots = task
    a = np.zeros(DIM, dtype=np.int64)
    b = a.copy()
    a[flat[~officer]] = (1 - 2 * x)[~officer]
    b[flat[officer]] = (1 - 2 * x)[officer]
    ja, jb = jac(a), jac(b)
    ta, tb = tangent(a), tangent(b)
    require(np.array_equal(quad(a) + 2 * quad(b), TARGET), 'base AME constant coefficient')
    require(not np.any(quad(a + b) - quad(a) - quad(b)), 'base AME sqrt2 coefficient')
    require(not np.any(ja @ ta + 2 * jb @ tb), 'JT constant coefficient')
    require(not np.any(ja @ tb + jb @ ta), 'JT sqrt2 coefficient')
    ws = []
    for small in rational_kernel(np.vstack([ja[:, other_cols], jb[:, other_cols]])):
        w = np.zeros(DIM, dtype=np.int64)
        w[other_cols] = small
        # With these bounds, every integer dot product below is far inside
        # int64 (36 terms of size at most 16^2 in each quadratic entry).
        require(np.max(np.abs(w)) <= 16, 'kernel coefficients need a fresh overflow bound')
        require(not np.any(ja @ w) and not np.any(jb @ w), 'extra kernel identity')
        ws.append(w)
    k = len(ws)
    require(k in (0, 2), f'unsupported candidate kernel dimension {k}')
    if k == 2:
        require(not np.any(quad(ws[0] + ws[1]) - quad(ws[0]) - quad(ws[1])),
                'nonzero mixed quadratic term: needs separate analysis')
    sa, sb = tangent(a, True), tangent(b, True)
    records = []
    for p, root in zip(PRIMES, roots):
        reduce_pair = lambda aa, bb: (aa % p + root * (bb % p) % p) % p
        j, t = reduce_pair(ja, jb), reduce_pair(ta, tb)
        rt, rj = rank_mod_p(t, p), rank_mod_p(j, p)
        require(rt == 60, f'orbit tangent rank {rt}')
        require(rj == DIM - 60 - k, f'Jacobian rank {rj} inconsistent with kernel ansatz k={k}')
        record = {'p': p, 'sqrt2': root, 'rank_T': rt, 'rank_J': rj}
        if k:
            rw = rank_mod_p(np.column_stack([t] + ws), p)
            rq = rank_mod_p(np.column_stack([j] + [quad(w) % p for w in ws]), p)
            require(rw == 62 and rq == rj + 2, 'second-order obstruction rank')
            record.update(rank_T_w=rw, rank_J_Q=rq)
        rs = rank_mod_p(reduce_pair(sa, sb), p)
        require(rs == 81, f'symmetric local action rank {rs}')
        record['rank_symmetric_local_action'] = rs
        records.append(record)
    return {'face': tag, 'diagonal_LU_orbit': oi, 'real_gauge_coset': coset,
            'signs': ''.join(map(str, x)), 'k': k,
            'kernel_vectors': [{str(i): int(w[i]) for i in np.flatnonzero(w)} for w in ws],
            'mixed_term_zero': True if k else None, 'primes': records,
            'certified_isolated': True, 'effective_stabiliser_lie_dim': 0}


def candidates(roots):
    tasks, coverage = [], []
    for tag in '+-':
        f = rp.Face(tag)
        other = rp.Face('-' if tag == '+' else '+')
        signs, solution_dim, gauge_dim = f.real_points()
        perms = f.automorphisms()
        orbits = f.orbits(signs, perms)
        g = fam.gauge_matrix(f.qs)
        kernel = gf2_nullspace(f.Pc)
        sign_rank = gf2_rank(g)
        basis = list(g.T)
        extra = []
        rank = sign_rank
        for v in kernel:
            new_rank = gf2_rank(np.array(basis + [v]))
            if new_rank > rank:
                extra.append(v)
                basis.append(v)
                rank = new_rank
        require(solution_dim == 32 and gauge_dim == 21 and len(kernel) == 21, 'gauge dimensions')
        require(not np.any(f.Pc @ g % 2), 'local signs belong to real gauge kernel')
        require(not np.any(f.Mg @ kernel.T % 2), 'real gauge preserves AME congruences')
        require(rank == 21 and len(extra) == (3 if tag == '+' else 4), 'full real gauge decomposition')
        require(sorted(i for orb in orbits for i in orb) == list(range(2048)), 'quotient orbit coverage')
        require(len(perms) == 192 and len(orbits) == (16 if tag == '+' else 24), 'support automorphisms')
        for img in perms:
            require(np.array_equal(np.sort(img), np.arange(60)), 'support permutation')
            require(np.array_equal(f.officer[img], f.officer), 'officer preservation')
            transformed = np.zeros_like(kernel)
            transformed[:, img] = kernel
            require(not np.any(f.Pc @ transformed.T % 2), 'gauge kernel covariant under support maps')
        offsets = [sum((bit * v for bit, v in zip(bits, extra)), np.zeros(60, dtype=np.int64)) % 2
                   for bits in itertools.product((0, 1), repeat=len(extra))]
        for oi, orbit in enumerate(orbits):
            for ci, offset in enumerate(offsets):
                x = (signs[orbit[0]] + offset) % 2
                require(not np.any((f.Mg @ x - 1) % 2), 'expanded candidate solves sign equations')
                tasks.append((tag, oi, ci, x, np.array(f.FLAT), f.officer,
                              np.array(other.FLAT)[other.officer], roots))
        coverage.append({'face': tag, 'sign_solution_dimension': solution_dim,
                         'real_diagonal_gauge_dimension': gauge_dim, 'local_sign_rank': sign_rank,
                         'extra_gauge_basis': [''.join(map(str, v)) for v in extra],
                         'extra_gauge_cosets': len(offsets), 'diagonal_LU_points': len(signs),
                         'diagonal_LU_orbits': len(orbits), 'candidates': len(orbits) * len(offsets)})
    return tasks, coverage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=min(4, os.cpu_count() or 1))
    args = parser.parse_args()
    if not 1 <= args.workers <= 15:
        parser.error('--workers must be in 1..15 (parent plus workers <=16)')
    start = time.monotonic()
    roots = []
    for p in PRIMES:
        require(p % 8 == 7 and all(p % d for d in range(2, int(p**0.5) + 1)), 'prime verification')
        root = pow(2, (p + 1) // 4, p)
        require(root * root % p == 2, 'finite-field sqrt2')
        roots.append(root)
    tasks, coverage = candidates(roots)
    for record in coverage:
        print('coverage ' + json.dumps(record, sort_keys=True), flush=True)
    check('orthogonal coverage: full real gauge = local signs + enumerated cosets', len(tasks) == 512)
    output = Path(__file__).resolve().parents[1] / 'results' / '16_real_isolation.jsonl'
    output.parent.mkdir(exist_ok=True)
    records = []
    try:
        with output.open('w') as stream:
            for record in coverage:
                stream.write(json.dumps({'coverage': record}, sort_keys=True) + '\n')
            if args.workers == 1:
                iterator = map(certify, tasks)
                pool = None
            else:
                # fork inherits read-only numpy/module state and avoids spawning extra
                # tracker processes; each task's arithmetic is single-threaded.
                pool = ProcessPoolExecutor(max_workers=args.workers,
                                           mp_context=multiprocessing.get_context('fork'))
                iterator = pool.map(certify, tasks, chunksize=1)
            try:
                for task, record in zip(tasks, iterator):
                    records.append(record)
                    stream.write(json.dumps(record, sort_keys=True) + '\n')
                    stream.flush()
                    print(f"[PASS] face {task[0]} orbit {task[1]:02d} gauge coset {task[2]:02d}: "
                          f"k={record['k']}, exact isolation; symmetric-action rank 81", flush=True)
            finally:
                if pool is not None:
                    pool.shutdown(wait=True, cancel_futures=True)
    except (ValueError, AssertionError) as exc:
        # Keep partial records and refuse a global certificate on an unexpected case.
        check('expanded certificate (stop on unsupported candidate)', False,
              f'candidate {len(records)}: {exc}')
        print(f'Runtime: {time.monotonic() - start:.1f} s', flush=True)
        done()
    for tag in '+-':
        rs = [r for r in records if r['face'] == tag]
        print(f"face {tag}: {len(rs)} candidates; k=0: {sum(r['k'] == 0 for r in rs)}; "
              f"k=2: {sum(r['k'] == 2 for r in rs)}", flush=True)
    check('all 512 expanded candidates certified isolated over Z[sqrt2]', len(records) == 512)
    check('effective projective local stabiliser Lie algebra zero (exact)',
          len(records) == 512 and all(r['effective_stabiliser_lie_dim'] == 0 for r in records))
    print(f'Runtime: {time.monotonic() - start:.1f} s ({args.workers} workers)', flush=True)
    done()


if __name__ == '__main__':
    main()
