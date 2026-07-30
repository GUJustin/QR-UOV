#!/usr/bin/env python3
"""Companion audit for projective target-orthogonal slicing and the QR-UOV scissors.

This script is intentionally self-contained except for importing the companion
`QRUOV_directional_audit.py` placed beside it. It:
  * checks the Chevalley--Warning feasibility inequalities;
  * computes fixed-oil seed-subspace and nonsingular-root success bounds;
  * projectivizes the nonzero-target residual and recomputes both cost ledgers;
  * certifies every sparse auxiliary-field polynomial used in the paper by
    Rabin's irreducibility criterion; and
  * computes attack-specific parameter thresholds for the combined scissors.

The cost rows are symbolic-algorithm stress estimates, not benchmark results.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "QRUOV_directional_audit.py"
spec = importlib.util.spec_from_file_location("directional_audit", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot load {BASE_PATH}")
base = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = base
spec.loader.exec_module(base)

Q = 127
LOGQ = math.log2(Q)
ELL = 3
PARAMS = [
    ("I", 156, 54, 52, 18, 143),
    ("III", 228, 78, 76, 26, 207),
    ("V", 306, 105, 102, 35, 272),
]

# All coefficients are +/-1. The map is exponent -> coefficient below X^d.
EXTRA_MODULI: Dict[int, Dict[int, int]] = {
    15: {12: 1, 1: -1, 0: 1},     # projective direct Level I
    16: {1: 1, 0: -1},             # direct Level I
    17: {8: 1, 0: 1},              # Level I +16 direct threshold
    22: {1: 1, 0: 1},              # direct Level III
    23: {2: 1, 0: 1},              # Level III direct target threshold
    25: {3: -1, 0: 1},             # intermediate Level III threshold
    26: {3: 1, 1: -1, 0: 1},       # Level III +16 direct threshold
    32: {2: 1, 0: -1},             # Level V direct target threshold
    34: {16: 1, 0: 1},             # Level V +16 direct threshold
}
base.SPARSE_MODULI.update(EXTRA_MODULI)


def log2sum(a: float, b: float) -> float:
    if a < b:
        a, b = b, a
    if b == -math.inf:
        return a
    return a + math.log2(1.0 + 2.0 ** (b - a))


def p_nonsingular_lower(r: int, q: int = Q) -> float:
    """Second-moment lower bound with dim(U cap oil) <= r-1."""
    mu = 1.0 - 1.0 / q
    pi = 1.0
    for j in range(1, r):
        pi *= 1.0 - q ** (-j)
    return mu * pi * pi / (2.0 + mu - (q - 1) * q ** (-r))


def log2_qpow_minus_one(a: int, q: int = Q) -> float:
    """Return log2(q^a-1) without materializing q^a."""
    tiny = q ** (-a) if a < 150 else 0.0
    return a * math.log2(q) + math.log2(1.0 - tiny)


def seed_intersection_failure_log2(v: int, m: int, seed_dim: int = 7, q: int = Q) -> float:
    """Union-bound log2 probability that a random seed subspace meets oil.

    For a uniform seed_dim-dimensional subspace S of F_q^(v+m) and any fixed
    m-dimensional oil subspace O,

      Pr[S cap O != {0}] <= |P(O)| * Pr[a fixed projective point lies in S]
                         = ((q^m-1)/(q-1))*((q^seed_dim-1)/(q^(v+m)-1)).

    This bound is over public attack coins and is valid for every fixed oil
    space.  It avoids converting an average anchor density into an unjustified
    expected-runtime statement for each fixed key.
    """
    n = v + m
    return (log2_qpow_minus_one(m, q) - math.log2(q - 1.0)
            + log2_qpow_minus_one(seed_dim, q)
            - log2_qpow_minus_one(n, q))


def seed_disjoint_probability_lower(v: int, m: int, seed_dim: int = 7, q: int = Q) -> float:
    log_fail = seed_intersection_failure_log2(v, m, seed_dim, q)
    if log_fail < -1000.0:
        return 1.0
    return 1.0 - 2.0 ** log_fail

def min_direct_extension_degree(k: int, q: int = Q) -> int:
    # |k'| >= 8*(2^k)^2 for the k-variable dehomogenized projective core.
    return math.ceil((3.0 + 2.0 * k) / math.log2(q) - 1e-15)


def sparse_extension_gate_envelope(d: int) -> float:
    return base.sparse_extension_gate_envelope(d)


def invertibility_probability_from_log(n: int, qlog: float) -> float:
    out = 1.0
    for j in range(1, n + 1):
        out *= 1.0 - 2.0 ** (-j * qlog)
    return out


def preprocessing_log2(v: int, m: int, r: int, restart: float) -> float:
    """Loose Boolean-gate envelope for seed search plus CW construction."""
    n = v + m
    gmul = 2.0 * LOGQ * LOGQ + LOGQ
    gadd = 4.0 * LOGQ

    projective_points = (Q ** 7 - 1) // (Q - 1)
    # One seed search and r-1 extension searches: r searches in total.
    # Each evaluates three seven-variable quadrics at every projective point.
    per_point = 3.0 * 7.0 * 7.0 * (gmul + gadd)
    cw = max(1, r) * projective_points * per_point

    # Generate a random full-rank seed basis and recompute all orthogonal
    # spaces by dense base-field elimination.
    linear = 16.0 * max(1, r) * n ** 3 * (gmul + gadd)
    return math.log2(restart * (cw + linear))


@dataclass
class DirectRow:
    level: str
    v: int
    m: int
    isotropic_dimension: int
    projective_core: int
    target: int
    cw_rhs: int
    cw_slack: int
    working_degree: int
    sparse_modulus: str
    seed_dimension: int
    seed_oil_intersection_failure_upper_log2: float
    seed_disjoint_probability_lower: float
    nonsingular_root_probability: float
    chart_probability: float
    separator_probability: float
    restart_factor: float
    lifting_stress: float
    auxiliary_stress: float
    preprocessing_stress: float
    one_flat_array_memory_log2_eib: float
    total_stress: float
    margin: float


def direct_row(level: str, v: int, m: int, V: int, O: int, target: int) -> DirectRow:
    del V, O
    n = v + m
    r = m - ELL
    k = r - 1
    rhs = (ELL + 1) * r + ELL
    if n < rhs:
        raise AssertionError("official profile fails Chevalley--Warning inequality")
    d = min_direct_extension_degree(k)
    if d not in base.SPARSE_MODULI:
        raise KeyError(f"no audited sparse modulus of degree {d}")
    qlog = d * LOGQ
    p_sep = 1.0 - 2.0 ** (2.0 * k - qlog)
    if p_sep <= 0:
        raise AssertionError("separator field too small")
    p_root = p_nonsingular_lower(r)
    p_chart = 1.0 - 1.0 / Q
    log_seed_fail = seed_intersection_failure_log2(v, m)
    p_seed = seed_disjoint_probability_lower(v, m)
    restart = 1.0 / (p_seed * p_root * p_chart * p_sep)

    log_c = float(k)
    log_c0 = k - 1.0 + math.log2(k + 2.0)
    log_outer = math.log2(k ** 3 + k ** 2)
    log_schedule_n = 3.0 + log_c + log_c0
    lifting = (log_outer + math.log2(sparse_extension_gate_envelope(d))
               + math.log2(restart) + math.log2(3.0)
               + base.log2_m_from_log2_n(log_schedule_n))

    # Explicit lower-stage envelope; the first term is quadratic factorization,
    # the second overcharges filtering and verification for every geometric root.
    # The projective resolution has at most 2^k points; each eligible point
    # gives at most two affine scalings, so the all-candidate verification term
    # remains bounded by 2^(k+1)=2^r.
    log_aux_inner = log2sum(2.0 * k,
                            math.log2(16.0 * m * n ** 3) + r)
    auxiliary = log_aux_inner + math.log2(sparse_extension_gate_envelope(d))
    prep = preprocessing_log2(v, m, r, restart)
    log_mem_bytes = log_c + 1.0 + log_c0 + math.log2(d * LOGQ) - 3.0
    total = log2sum(log2sum(lifting, auxiliary), prep)

    return DirectRow(
        level=level, v=v, m=m, isotropic_dimension=r,
        projective_core=k, target=target,
        cw_rhs=rhs, cw_slack=n-rhs, working_degree=d,
        sparse_modulus=base.modulus_string(d),
        seed_dimension=7,
        seed_oil_intersection_failure_upper_log2=log_seed_fail,
        seed_disjoint_probability_lower=p_seed,
        nonsingular_root_probability=p_root,
        chart_probability=p_chart,
        separator_probability=p_sep, restart_factor=restart,
        lifting_stress=lifting, auxiliary_stress=auxiliary,
        preprocessing_stress=prep,
        one_flat_array_memory_log2_eib=log_mem_bytes - 60.0,
        total_stress=total,
        margin=target-total,
    )


def find_direct_threshold(start_m: int, target_cost: float) -> tuple[int, float, int]:
    for m in range(start_m, 300):
        if m % ELL != 0:
            continue
        r = m - ELL
        k = r - 1
        d = min_direct_extension_degree(k)
        if d not in base.SPARSE_MODULI:
            continue
        # The v value affects only a negligible seed-intersection term; use a generous v.
        row = direct_row("threshold", 3 * max(r, 1), m, max(r, 1), 1, 0)
        if row.total_stress >= target_cost:
            return m, row.total_stress, d
    raise RuntimeError("direct threshold not found")


def find_directional_threshold(start_V: int, target_cost: float) -> tuple[int, float, int]:
    for V in range(start_V, 300):
        d = 3 * base.min_extension_degree(V)
        if d not in base.SPARSE_MODULI:
            continue
        row = base.cost_row("threshold", V, 1, 1, 0)
        if row.corrected_sparse_field_stress >= target_cost:
            return V, row.corrected_sparse_field_stress, d
    raise RuntimeError("directional threshold not found")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", type=Path, default=Path("QRUOV_scissors_audit.json"))
    args = ap.parse_args()

    direct = [direct_row(*p) for p in PARAMS]
    directional = [base.cost_row(level, V, O, m, target)
                   for level, v, m, V, O, target in PARAMS]

    cert_degrees = sorted({r.working_degree for r in direct}
                          | {r.d for r in directional}
                          | {17, 23, 25, 26, 27, 32, 33, 34, 36})
    certs = [base.rabin_irreducibility_certificate(d) for d in cert_degrees]

    combined = []
    for p, dr, kr in zip(PARAMS, direct, directional):
        level, v, m, V, O, target = p
        best = min(dr.total_stress, kr.corrected_sparse_field_stress)
        source = "target-orthogonal direct forgery" if dr.total_stress <= kr.corrected_sparse_field_stress else "directional equivalent-key recovery"
        meet_V0 = find_directional_threshold(V, target)
        meet_m = find_direct_threshold(m, target)
        meet_pair_V = max(meet_V0[0], meet_m[0] - ELL)
        meet_V_row = base.cost_row("threshold", meet_pair_V, 1, 1, 0)
        meet_V = (meet_pair_V, meet_V_row.corrected_sparse_field_stress, meet_V_row.d)

        cushion_V0 = find_directional_threshold(V, target + 16)
        cushion_m = find_direct_threshold(m, target + 16)
        cushion_pair_V = max(cushion_V0[0], cushion_m[0] - ELL)
        cushion_V_row = base.cost_row("threshold", cushion_pair_V, 1, 1, 0)
        cushion_V = (cushion_pair_V, cushion_V_row.corrected_sparse_field_stress, cushion_V_row.d)

        if meet_m[0] % ELL or cushion_m[0] % ELL:
            raise AssertionError("QR-UOV output count must be divisible by quotient degree")
        if meet_V[0] < meet_m[0] - ELL or cushion_V[0] < cushion_m[0] - ELL:
            raise AssertionError("reported pair does not support the s=ell direct reduction")

        combined.append({
            "level": level,
            "current_V": V,
            "current_m": m,
            "target": target,
            "direct_cost": dr.total_stress,
            "directional_cost": kr.corrected_sparse_field_stress,
            "best_cost": best,
            "best_attack": source,
            "best_margin": target - best,
            "simple_qr_compatible_target_pair": {
                "V": meet_V[0], "directional_cost": meet_V[1], "directional_field_degree": meet_V[2],
                "m": meet_m[0], "direct_cost": meet_m[1], "direct_field_degree": meet_m[2],
            },
            "simple_qr_compatible_target_plus_16_pair": {
                "V": cushion_V[0], "directional_cost": cushion_V[1], "directional_field_degree": cushion_V[2],
                "m": cushion_m[0], "direct_cost": cushion_m[1], "direct_field_degree": cushion_m[2],
            },
        })

    # Frozen headline checks.
    expected_direct = {"I": 151.639, "III": 203.507, "V": 260.389}
    for row in direct:
        if abs(row.total_stress - expected_direct[row.level]) > 5e-4:
            raise AssertionError(f"direct headline drift at Level {row.level}: {row.total_stress}")
    expected_best = {"I": 151.639, "III": 203.507, "V": 260.389}
    for row in combined:
        if abs(row["best_cost"] - expected_best[row["level"]]) > 5e-4:
            raise AssertionError(f"combined headline drift at Level {row['level']}")

    payload = {
        "description": "Exact projective slicing inequalities, fixed-oil seed bounds, Karatsuba sparse-field certificates, and combined QR-UOV attack ledger.",
        "direct_forgery_rows": [asdict(r) for r in direct],
        "directional_recovery_rows": [asdict(r) for r in directional],
        "combined_scissors_rows": combined,
        "sparse_modulus_certificates": certs,
        "scissors_identity": {
            "quotient_degree": ELL,
            "statement": "For nonzero targets at fixed m, directional recovery covers V <= m-ell-1 and projectivized target-orthogonal slicing covers V >= m-ell, so one attack core is at most m-ell-1.",
        },
    }
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print("QR-UOV target-orthogonal slicing and scissors audit")
    print("=" * 72)
    for dr, kr, row in zip(direct, directional, combined):
        print(
            f"Level {dr.level}: CW slack={dr.cw_slack}, isotropic dim={dr.isotropic_dimension}, "
            f"projective core={dr.projective_core}, "
            f"log2(seed failure)<={dr.seed_oil_intersection_failure_upper_log2:.3f}, "
            f"p_root>={dr.nonsingular_root_probability:.12f}, d={dr.working_degree}, "
            f"direct={dr.total_stress:.3f}, directional={kr.corrected_sparse_field_stress:.3f}, "
            f"best={row['best_cost']:.3f} ({row['best_margin']:+.3f} target margin)"
        )
        t = row["simple_qr_compatible_target_pair"]
        c = row["simple_qr_compatible_target_plus_16_pair"]
        print(
            f"  target pair: V={t['V']} (cost {t['directional_cost']:.3f}), "
            f"m={t['m']} (cost {t['direct_cost']:.3f}); "
            f"target+16 pair: V={c['V']} (cost {c['directional_cost']:.3f}), "
            f"m={c['m']} (cost {c['direct_cost']:.3f})"
        )
    for c in certs:
        print("Rabin certificate:", c["polynomial"], c["gcd_degrees_for_prime_divisors"])
    print("Wrote", args.json)


if __name__ == "__main__":
    main()
