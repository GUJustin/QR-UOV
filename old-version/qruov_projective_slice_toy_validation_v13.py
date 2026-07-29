#!/usr/bin/env python3
"""Exact toy validation of random projective slicing for central UOV forms.

This is not evidence at production dimensions.  It checks, over small prime
fields, the new end-to-end identities used by the random-projective-slice
attack:
  * a public V-plane meets the planted oil space in a linear point;
  * a smooth full core has that point as a regular root;
  * the ambient derivative kernel at that root is exactly the oil space.

For V=2 the smoothness test is exact over the algebraic closure: on each
projective chart we compute a Groebner basis of the two restricted conics and
their affine Jacobian determinant.  A unit Groebner basis on every chart is
equivalent to absence of geometric singular roots.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from typing import Iterable

import sympy as sp


Matrix = list[list[int]]


def matmul(a: Matrix, b: Matrix, p: int) -> Matrix:
    assert len(a[0]) == len(b)
    return [[sum(a[i][k] * b[k][j] for k in range(len(b))) % p
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a: Matrix) -> Matrix:
    return [list(row) for row in zip(*a)]


def matadd(a: Matrix, b: Matrix, p: int) -> Matrix:
    return [[(x + y) % p for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def matscale(c: int, a: Matrix, p: int) -> Matrix:
    return [[c * x % p for x in row] for row in a]


def rank_mod(a: Matrix, p: int) -> int:
    if not a:
        return 0
    m = [row[:] for row in a]
    rows, cols = len(m), len(m[0])
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if m[i][c] % p), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        inv = pow(m[r][c] % p, -1, p)
        m[r] = [(inv * x) % p for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] % p:
                f = m[i][c] % p
                m[i] = [(x - f * y) % p for x, y in zip(m[i], m[r])]
        r += 1
        if r == rows:
            break
    return r


def nullspace_mod(a: Matrix, p: int) -> list[list[int]]:
    """Return a basis of the right nullspace of a over F_p."""
    if not a:
        return []
    m = [row[:] for row in a]
    rows, cols = len(m), len(m[0])
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if m[i][c] % p), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        inv = pow(m[r][c] % p, -1, p)
        m[r] = [(inv * x) % p for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] % p:
                f = m[i][c] % p
                m[i] = [(x - f * y) % p for x, y in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis: list[list[int]] = []
    for f in free:
        v = [0] * cols
        v[f] = 1
        for i, c in enumerate(pivots):
            v[c] = (-m[i][f]) % p
        basis.append(v)
    return basis


def same_span(columns_a: list[list[int]], columns_b: list[list[int]], p: int) -> bool:
    if not columns_a or not columns_b:
        return columns_a == columns_b
    a = transpose(columns_a)  # columns supplied as vectors -> matrix columns
    b = transpose(columns_b)
    ra = rank_mod(a, p)
    rb = rank_mod(b, p)
    joined = [ra_row + rb_row for ra_row, rb_row in zip(a, b)]
    return ra == rb == rank_mod(joined, p)


def random_matrix(rows: int, cols: int, p: int, rng: random.Random) -> Matrix:
    return [[rng.randrange(p) for _ in range(cols)] for _ in range(rows)]


def random_symmetric(n: int, p: int, rng: random.Random) -> Matrix:
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            x = rng.randrange(p)
            a[i][j] = a[j][i] = x
    return a


def central_matrix(a: Matrix, b: Matrix, p: int) -> Matrix:
    v, o = len(a), len(b[0])
    out = [[0] * (v + o) for _ in range(v + o)]
    for i in range(v):
        for j in range(v):
            out[i][j] = a[i][j] % p
    for i in range(v):
        for j in range(o):
            out[i][v + j] = b[i][j] % p
            out[v + j][i] = b[i][j] % p
    return out


def quadratic_expr(h: Matrix, vars_: list[sp.Symbol], p: int) -> sp.Expr:
    n = len(vars_)
    expr = 0
    for i in range(n):
        expr += h[i][i] * vars_[i] ** 2
        for j in range(i + 1, n):
            expr += 2 * h[i][j] * vars_[i] * vars_[j]
    return sp.Poly(expr, *vars_, modulus=p).as_expr()


def smooth_two_conic_core(hs: list[Matrix], p: int) -> bool:
    """Exact geometric smoothness test for two conics in P^2 over F_p."""
    assert len(hs) == 2 and len(hs[0]) == 3
    t = sp.symbols("t0:3")
    polys = [quadratic_expr(h, list(t), p) for h in hs]
    for chart in range(3):
        free_idx = [i for i in range(3) if i != chart]
        free = [t[i] for i in free_idx]
        subs = {t[chart]: 1}
        f1 = sp.Poly(polys[0].subs(subs), *free, modulus=p).as_expr()
        f2 = sp.Poly(polys[1].subs(subs), *free, modulus=p).as_expr()
        det = sp.diff(f1, free[0]) * sp.diff(f2, free[1]) - sp.diff(f1, free[1]) * sp.diff(f2, free[0])
        det = sp.Poly(det, *free, modulus=p).as_expr()
        gb = sp.groebner([f1, f2, det], *free, modulus=p)
        if not any(poly.as_expr() == 1 for poly in gb.polys):
            return False
    return True


def combine_matrices(coeffs: list[int], mats: list[Matrix], p: int) -> Matrix:
    n = len(mats[0])
    out = [[0] * n for _ in range(n)]
    for c, a in zip(coeffs, mats):
        out = matadd(out, matscale(c, a, p), p)
    return out


@dataclass
class Stats:
    keys: int = 0
    slice_trials: int = 0
    smooth_cores: int = 0
    recovered: int = 0
    rank_deficient_intersections: int = 0


def one_key(V: int, O: int, m: int, p: int, slice_trials: int, rng: random.Random, stats: Stats) -> None:
    assert V == 2, "The exact Groebner smoothness checker currently supports V=2."
    n = V + O
    forms: list[Matrix] = []
    for _ in range(m):
        forms.append(central_matrix(random_symmetric(V, p, rng), random_matrix(V, O, p, rng), p))
    stats.keys += 1

    oil_basis = [[0] * V + [1 if j == a else 0 for j in range(O)] for a in range(O)]

    for _ in range(slice_trials):
        stats.slice_trials += 1
        # Public projective V-plane: an N x (V+1) full-column-rank matrix.
        for _attempt in range(100):
            M = random_matrix(n, V + 1, p, rng)
            if rank_mod(M, p) == V + 1:
                break
        else:
            raise RuntimeError("failed to sample full-rank slice")

        # The plane-oil intersection solves top(M) t = 0.
        top = M[:V]
        ns = nullspace_mod(top, p)
        if len(ns) != 1:
            stats.rank_deficient_intersections += 1
            continue
        tstar = ns[0]
        xstar = [sum(M[i][j] * tstar[j] for j in range(V + 1)) % p for i in range(n)]
        if any(xstar[i] % p for i in range(V)) or not any(xstar[V + a] % p for a in range(O)):
            raise AssertionError("computed intersection is not a nonzero oil point")

        R = random_matrix(V, m, p, rng)
        combined = [combine_matrices(R[a], forms, p) for a in range(V)]
        restricted = [matmul(transpose(M), matmul(h, M, p), p) for h in combined]
        if not smooth_two_conic_core(restricted, p):
            continue
        stats.smooth_cores += 1

        # Every central form vanishes at xstar.
        for h in forms:
            hx = [sum(h[i][j] * xstar[j] for j in range(n)) % p for i in range(n)]
            val = sum(xstar[i] * hx[i] for i in range(n)) % p
            if val:
                raise AssertionError("oil intersection does not satisfy a central form")

        # C_x has columns P_i x; its transpose kernel should be exactly oil.
        ccols = []
        for h in forms:
            ccols.append([sum(h[i][j] * xstar[j] for j in range(n)) % p for i in range(n)])
        C = transpose(ccols)  # N x m
        ker = nullspace_mod(transpose(C), p)
        if len(ker) != O or not same_span(ker, oil_basis, p):
            raise AssertionError("derivative kernel failed to recover the oil space")
        stats.recovered += 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", type=int, default=7, help="odd prime field size")
    ap.add_argument("--keys", type=int, default=20)
    ap.add_argument("--trials", type=int, default=10, help="slice/output trials per key")
    ap.add_argument("--oil", type=int, default=2)
    ap.add_argument("--equations", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20260722)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument(
        "--suite",
        dest="suite",
        action="store_true",
        default=True,
        help="run the frozen three-shape audit suite (default)",
    )
    mode.add_argument(
        "--single",
        dest="suite",
        action="store_false",
        help="run one custom case using --field/--oil/--equations/--keys/--trials/--seed",
    )
    args = ap.parse_args()

    cases = [(args.field, args.oil, args.equations, args.keys, args.trials, args.seed)]
    if args.suite:
        cases = [
            (5, 2, 4, 20, 12, 1),
            (7, 3, 4, 15, 10, 2),
            (11, 2, 3, 12, 10, 3),
        ]

    total = Stats()
    for field, oil, equations, keys, trials, seed in cases:
        if field % 2 == 0:
            raise SystemExit("field must have odd characteristic")
        rng = random.Random(seed)
        stats = Stats()
        for _ in range(keys):
            one_key(2, oil, equations, field, trials, rng, stats)
        print("Exact V=2 projective-slice validation")
        print(f"field=F_{field}, O={oil}, m={equations}, keys={stats.keys}")
        print(f"slice/output trials: {stats.slice_trials}")
        print(f"rank-deficient plane-oil intersections: {stats.rank_deficient_intersections}")
        print(f"geometrically smooth full cores: {stats.smooth_cores}")
        print(f"successful derivative-kernel recoveries: {stats.recovered}")
        if stats.smooth_cores != stats.recovered:
            raise SystemExit("FAIL: a smooth core did not recover the oil space")
        print("PASS: every smooth core recovered the exact planted oil space\n")
        for name in total.__dataclass_fields__:
            setattr(total, name, getattr(total, name) + getattr(stats, name))

    if args.suite:
        print("Frozen suite aggregate")
        print(f"keys: {total.keys}")
        print(f"slice/output trials: {total.slice_trials}")
        print(f"geometrically smooth full cores: {total.smooth_cores}")
        print(f"successful derivative-kernel recoveries: {total.recovered}")
        print("PASS: all suite checks succeeded")


if __name__ == "__main__":
    main()
