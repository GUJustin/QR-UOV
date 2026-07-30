#!/usr/bin/env python3
"""Independent recomputation of the QR-UOV split-before-lift headline ledger.

This file deliberately imports neither the primary split-lift audit nor the parent
audit modules.  It duplicates the formulas from the technical note in a different
code structure and checks the published JSON headline values.
"""
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
    15: {12: 1, 1: -1, 0: 1},
    22: {1: 1, 0: 1},
    30: {7: 1, 0: 1},
}
ROWS = [
    ("I", 156, 54, 143),
    ("III", 228, 78, 207),
    ("V", 306, 105, 272),
]


def mulcost(n: int) -> float:
    return float(n * n) if n < 4 else n * math.log2(n) * math.log2(math.log2(n))


def lsum(vals: list[float]) -> float:
    m = max(vals)
    return m + math.log2(sum(2 ** (v - m) for v in vals))


def log_qpow_minus_one(e: int) -> float:
    tiny = Q ** (-e) if e < 150 else 0.0
    return e * LQ + math.log2(1.0 - tiny)


def seed_success(v: int, m: int) -> float:
    n = v + m
    lf = (log_qpow_minus_one(m) - math.log2(Q - 1)
          + log_qpow_minus_one(7) - log_qpow_minus_one(n))
    return 1.0 - 2 ** lf


@functools.lru_cache(None)
def conv_plan(d: int) -> tuple[float, int, int]:
    sm, sa = d * d, (d - 1) ** 2
    best = (sm * GMUL + sa * GADD, sm, sa)
    if d >= 2:
        hi, lo = (d + 1) // 2, d // 2
        ph, pl = conv_plan(hi), conv_plan(lo)
        km = 2 * ph[1] + pl[1]
        ka = 2 * ph[2] + pl[2] + 4 * d - 4
        cand = (km * GMUL + ka * GADD, km, ka)
        if cand[0] < best[0]:
            best = cand
    return best


def extmul(d: int) -> float:
    return conv_plan(d)[0] + len(MODULI[d]) * (d - 1) * GADD


def relaxed(n: int) -> float:
    ell = math.ceil(math.log2(n))
    return (sum(2 * math.ceil(n / (1 << j)) * mulcost(1 << j)
                for j in range(ell + 1))
            + 8 * n * (ell + 1))


def long_ext_product(n: int, d: int) -> float:
    s = 2 * d - 1
    return (s * mulcost(n) * GMUL
            + 2 * n * s * d * (GMUL + GADD)
            + 2 * n * (s * s + 8 * d * d) * (GMUL + GADD))


def root_lower(r: int) -> float:
    mu = 1 - 1 / Q
    pi = math.prod(1 - Q ** (-j) for j in range(1, r))
    return mu * pi / 2 - mu * mu / 8 + (Q - 1) * mu * Q ** (-r) / 8


def preprocessing(v: int, m: int, r: int, restart: float) -> float:
    n = v + m
    points = (Q ** 7 - 1) // (Q - 1)
    cw = r * points * 3 * 7 * 7 * (GMUL + GADD)
    linear = 16 * r * n ** 3 * (GMUL + GADD)
    return restart * (cw + linear)


def compute(level: str, v: int, m: int, target: int) -> dict[str, float]:
    a, r = m - 4, m - 3
    d = math.ceil((3 + 2 * a) / LQ - 1e-15)
    C = 1 << a
    C0 = (1 << (a - 1)) * (a + 2)
    P = 4 * C0 + 1
    N = 4 * C * P
    B = 1 << 16

    failures = [C * C / Q ** d, 2 * C / Q ** d, C / Q ** d]
    p = seed_success(v, m) * (root_lower(r) - 1 / (2 * (B + 1))) * (1 - 1 / Q)
    for f in failures:
        p *= 1 - f
    restart = 1 / p

    h = a * (a + 1) // 2
    branch_one = (h * relaxed(P) * GMUL
                  + (a * h + a * a + 8 * a) * P * (GMUL + GADD))
    all_branches = C * branch_one
    branch = restart * all_branches
    filt_series = restart * C * 8 * mulcost(P) * extmul(d)

    tree_units = 7 * a / 2
    rr = (3 * C * 8 * mulcost(2 * C0) * math.log2(2 * C0)
          + C * 8 * mulcost(4 * C0) * math.log2(4 * C0)) / mulcost(N)
    candidates = 0.5 + B * C / Q ** d
    candidate_units = (1 + 2 * a) * C * mulcost(P) / mulcost(N)
    terminal_units = tree_units + rr + candidates * candidate_units + 64
    terminal = restart * terminal_units * long_ext_product(N, d)
    relifts = restart * candidates * all_branches
    leaves = restart * C * (3 * a + 16) * P * extmul(d)

    n = v + m
    aux_inner_log = lsum([2.0 * a, math.log2(16 * m * n ** 3) + r])
    aux = restart * (2 ** aux_inner_log) * extmul(d)
    prep = preprocessing(v, m, r, restart)

    logs = {
        "initial_base_field_branch_lifts": math.log2(branch),
        "valuation_safe_filter_series": math.log2(filt_series),
        "terminal_scalar_recombination": math.log2(terminal),
        "on_demand_coordinate_relifts": math.log2(relifts),
        "extension_leaf_formation": math.log2(leaves),
        "auxiliary_factorization_and_verification": math.log2(aux),
        "preprocessing": math.log2(prep),
    }
    logs["total"] = lsum(list(logs.values()))
    logs["margin"] = target - logs["total"]
    return logs


def main() -> None:
    here = Path(__file__).resolve().parent
    expected = json.loads((here / "QRUOV_split_lift_audit.json").read_text())
    exp_by_level = {r["level"]: r["cost_log2_gates"] for r in expected["rows"]}
    print("Independent split-before-lift ledger verification")
    print("=" * 60)
    for args in ROWS:
        level = args[0]
        got = compute(*args)
        exp = exp_by_level[level]
        for key in got:
            if abs(got[key] - exp[key]) > 5e-10:
                raise AssertionError((level, key, got[key], exp[key]))
        print(f"Level {level}: total={got['total']:.9f}, margin={got['margin']:+.9f} [MATCH]")


if __name__ == "__main__":
    main()
