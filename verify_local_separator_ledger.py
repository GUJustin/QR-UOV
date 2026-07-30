#!/usr/bin/env python3
"""Independent no-import verification of the local-separator v2 ledger."""
from __future__ import annotations

import functools
import json
import math
from pathlib import Path

Q = 127
LQ = math.log2(Q)
GMUL = 2 * LQ * LQ + LQ
GADD = 4 * LQ
MODULI = {
    8: {4: 1, 0: -1},
    9: {2: 1, 0: 1},
    12: {3: -1, 1: 1, 0: 1},
    16: {1: 1, 0: -1},
    15: {12: 1, 1: -1, 0: 1},
    22: {1: 1, 0: 1},
    30: {7: 1, 0: 1},
}
ROWS = [
    ("I", 156, 54, 143, 8, 15),
    ("III", 228, 78, 207, 12, 22),
    ("V", 306, 105, 272, 16, 30),
]


def M(n):
    return float(n * n) if n < 4 else n * math.log2(n) * math.log2(math.log2(n))


def lsum(vals):
    top = max(vals)
    return top + math.log2(sum(2 ** (v - top) for v in vals))


def log_qm1(e):
    tiny = Q ** (-e) if e < 150 else 0.0
    return e * LQ + math.log2(1 - tiny)


def seed_success(v, m):
    n = v + m
    lf = log_qm1(m) - math.log2(Q - 1) + log_qm1(7) - log_qm1(n)
    return 1 - 2**lf


@functools.lru_cache(None)
def conv(d):
    sm, sa = d * d, (d - 1) ** 2
    best = (sm * GMUL + sa * GADD, sm, sa)
    if d >= 2:
        hi, lo = (d + 1) // 2, d // 2
        ph, pl = conv(hi), conv(lo)
        km = 2 * ph[1] + pl[1]
        ka = 2 * ph[2] + pl[2] + 4 * d - 4
        cand = (km * GMUL + ka * GADD, km, ka)
        if cand[0] < best[0]:
            best = cand
    return best


def extmul(d):
    return conv(d)[0] + len(MODULI[d]) * (d - 1) * GADD


def relaxed(n):
    ell = math.ceil(math.log2(n))
    return sum(2 * math.ceil(n / (1 << j)) * M(1 << j) for j in range(ell + 1)) + 8 * n * (ell + 1)


def longprod(n, d):
    s = 2 * d - 1
    return (
        s * M(n) * GMUL
        + 2 * n * s * d * (GMUL + GADD)
        + 2 * n * (s * s + 8 * d * d) * (GMUL + GADD)
    )


def root_lower(r):
    mu = 1 - 1 / Q
    pi = math.prod(1 - Q ** (-j) for j in range(1, r))
    return mu * pi / 2 - mu * mu / 8 + (Q - 1) * mu * Q ** (-r) / 8


def preprocessing(v, m, r, restart):
    n = v + m
    points = (Q**7 - 1) // (Q - 1)
    return restart * (
        r * points * 3 * 7 * 7 * (GMUL + GADD)
        + 16 * r * n**3 * (GMUL + GADD)
    )


def compute(level, v, m, target, d, old_d):
    a, r = m - 4, m - 3
    C = 1 << a
    C0 = (1 << (a - 1)) * (a + 2)
    P = 4 * C0 + 1
    N = 4 * C * P
    K = Q**d
    p = seed_success(v, m) * root_lower(r) * (1 - 1 / Q) * (1 - (C - 1) / K) * (1 - C / (K - C))
    restart = 1 / p
    h = a * (a + 1) // 2
    one = h * relaxed(P) * GMUL + (a * h + a * a + 8 * a) * P * (GMUL + GADD)
    allb = C * one
    branch = restart * allb
    filt = restart * C * 8 * M(P) * extmul(d)
    tree = 5 * a / 2
    rr = (
        2 * C * 8 * M(2 * C0) * math.log2(2 * C0)
        + C * 8 * M(4 * C0) * math.log2(4 * C0)
    ) / M(N)
    cand = (1 + 2 * a) * C * M(P) / M(N)
    terminal = (restart * (tree + rr + 64) + cand) * longprod(N, d)
    relifts = allb
    leaves = restart * C * (2 * a + 16) * P * extmul(d)
    coordinate_leaves = C * a * P * extmul(d)
    n = v + m
    aux_inner = lsum([2 * a, math.log2(16 * m * n**3) + r])
    aux = restart * (2**aux_inner) * extmul(old_d)
    prep = preprocessing(v, m, r, restart)
    logs = {
        "initial_base_field_branch_lifts": math.log2(branch),
        "valuation_safe_filter_series": math.log2(filt),
        "terminal_scalar_recombination": math.log2(terminal),
        "on_demand_coordinate_relifts": math.log2(relifts),
        "extension_leaf_formation": math.log2(leaves),
        "coordinate_relift_leaf_formation": math.log2(coordinate_leaves),
        "auxiliary_factorization_and_verification": math.log2(aux),
        "preprocessing": math.log2(prep),
    }
    logs["total"] = lsum(list(logs.values()))
    logs["margin"] = target - logs["total"]
    return logs


def main():
    here = Path(__file__).resolve().parent
    data = json.loads((here / "QRUOV_local_separator_audit.json").read_text())
    expected = {r["level"]: r["cost_log2_gates"] for r in data["rows"]}
    print("Independent local-separator v2 ledger verification")
    print("=" * 64)
    for row in ROWS:
        got = compute(*row)
        want = expected[row[0]]
        for key, value in got.items():
            if abs(value - want[key]) > 5e-10:
                raise AssertionError((row[0], key, value, want[key]))
        print(f"Level {row[0]}: total={got['total']:.9f}, margin={got['margin']:+.9f} [MATCH]")


if __name__ == "__main__":
    main()
