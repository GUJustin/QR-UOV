#!/usr/bin/env python3
"""Retry-aware split-before-lift attack ledger for QR-UOV.

This audit keeps the Boolean-gate model and sparse-field certificates of the
strict QR-UOV bundle, but replaces global geometric-resolution lifting by:

  * independent base-field Hensel lifts of all Vandermonde start branches;
  * relaxed products for the quadratic feedback terms;
  * terminal scalar recombination over the separator field;
  * a two-scalar Frobenius descent test;
  * the valuation-safe oil filter H0/(c + lambda.x); and
  * on-demand coordinate recovery for the few surviving rational roots.

The output is a finite stress ledger, not a measured implementation.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
STRICT_AUDIT = HERE / "parent_audit" / "QRUOV_scissors_audit.py"
if not STRICT_AUDIT.exists():
    STRICT_AUDIT = Path(
        "/mnt/data/qruov_strict/QRUOV_strict_final_bundle/QRUOV_scissors_audit.py"
    )
spec = importlib.util.spec_from_file_location("strict_scissors", STRICT_AUDIT)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {STRICT_AUDIT}")
strict = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = strict
spec.loader.exec_module(strict)

Q = 127
LOGQ = math.log2(Q)
GMUL = 2.0 * LOGQ * LOGQ + LOGQ
GADD = 4.0 * LOGQ
RR_CONST = 8.0
FILTER_INVERSE_CONST = 8.0
MISC_LONG_PRODUCTS = 64.0

# Rounded, review-friendly caps.  Their exact optimization is immaterial.
ROOT_CAPS = {"I": 1 << 16, "III": 1 << 16, "V": 1 << 16}
PARAMS = [
    ("I", 156, 54, 52, 18, 143),
    ("III", 228, 78, 76, 26, 207),
    ("V", 306, 105, 102, 35, 272),
]


def M(n: int) -> float:
    """The same explicit multiplication stress function as the strict bundle."""
    if n < 4:
        return float(n * n)
    return n * math.log2(n) * math.log2(math.log2(n))


def log2sum(*logs: float) -> float:
    top = max(logs)
    return top + math.log2(sum(2.0 ** (x - top) for x in logs))


def log_success_from_failure(failure: float) -> float:
    if not 0.0 <= failure < 1.0:
        raise ValueError(f"invalid failure probability {failure}")
    return math.log1p(-failure)


def relaxed_upper(n: int) -> float:
    """Deliberately loose dyadic envelope for one online series product.

    Published online multiplication bounds are much smaller in the FFT model.
    This explicit sum is retained as a cushion and is charged in base-field
    multiplication units.
    """
    levels = math.ceil(math.log2(n))
    out = 0.0
    for j in range(levels + 1):
        block = 1 << j
        out += 2.0 * math.ceil(n / block) * M(block)
    out += 8.0 * n * (levels + 1)
    return out


def projective_root_lower(r: int) -> float:
    """Second-Bonferroni lower bound for an eligible nonsingular point."""
    mu = 1.0 - 1.0 / Q
    pi = math.prod(1.0 - Q ** (-j) for j in range(1, r))
    return mu * pi / 2.0 - mu * mu / 8.0 + (Q - 1) * mu * Q ** (-r) / 8.0


def long_extension_product_gates(n: int, d: int) -> float:
    """One long F_(127^d)[Z] product via 2d-1 base-field evaluations.

    Coefficients are degree-<d polynomials in the extension generator.  Evaluate
    at 2d-1 distinct F_127 points, perform 2d-1 base polynomial products,
    interpolate the degree-<2d-1 coefficient products, and reduce modulo the
    audited sparse irreducible polynomial.  The dense transforms below are
    intentionally overcharged.
    """
    points = 2 * d - 1
    if points > Q:
        raise ValueError("not enough base-field evaluation points")
    base_products = points * M(n) * GMUL
    evaluations = 2.0 * n * points * d * (GMUL + GADD)
    interpolation_and_reduction = (
        2.0 * n * (points * points + 8.0 * d * d) * (GMUL + GADD)
    )
    return base_products + evaluations + interpolation_and_reduction


@dataclass
class Row:
    level: str
    target: int
    a: int
    r: int
    extension_degree: int
    C: int
    C0: int
    precision: int
    root_cap: int
    probabilities: dict[str, float]
    multiplication: dict[str, float]
    terminal: dict[str, float]
    cost_log2_gates: dict[str, float]
    sensitivity: dict[str, float]


def attack_row(
    level: str,
    v_expanded: int,
    m: int,
    V: int,
    O: int,
    target: int,
    *,
    rr_const: float = RR_CONST,
    root_cap: int | None = None,
) -> Row:
    a = m - 4
    r = m - 3
    old = strict.direct_row(level, v_expanded, m, V, O, target)
    d = old.working_degree
    cap = ROOT_CAPS[level] if root_cap is None else root_cap

    C = 1 << a
    C0 = (1 << (a - 1)) * (a + 2)
    # Degree <=2C0 for the graph of H0/(c+lambda.x), plus one coefficient.
    precision = 4 * C0 + 1
    flat_length = 4 * C * precision
    quadratic_monomials = a * (a + 1) // 2

    p_root = projective_root_lower(r)
    p_root_cap = p_root - 1.0 / (2.0 * (cap + 1))
    p_chart = 1.0 - 1.0 / Q
    p_seed = strict.seed_disjoint_probability_lower(v_expanded, m)

    sep_failure = C * C / (Q**d)
    denominator_failure = 2.0 * C / (Q**d)
    frobenius_filter_failure = C / (Q**d)

    log_success = (
        math.log(p_seed)
        + math.log(p_root_cap)
        + math.log(p_chart)
        + log_success_from_failure(sep_failure)
        + log_success_from_failure(denominator_failure)
        + log_success_from_failure(frobenius_filter_failure)
    )
    restart = math.exp(-log_success)

    relaxed = relaxed_upper(precision)
    one_branch = (
        quadratic_monomials * relaxed * GMUL
        + (a * quadratic_monomials + a * a + 8 * a)
        * precision
        * (GMUL + GADD)
    )
    all_branches = C * one_branch
    branch_expected = restart * all_branches

    # h=H0/(c+lambda.x): eight M(P) extension operations per branch cover
    # denominator inversion and multiplication by H0.
    rational_filter_series = (
        restart
        * C
        * FILTER_INVERSE_CONST
        * M(precision)
        * strict.sparse_extension_gate_envelope(d)
    )

    # Three value families: rho.x, rho^(q^{-1}).x, and H0/(c+lambda.x).
    value_families = 3
    tree_units = (1.0 + 2.0 * value_families) * a / 2.0

    # q and two linear numerators have t-degree <=C0.  The rationalized H0
    # numerator has degree <=2C0.
    rr_linear = 3.0 * C * rr_const * M(2 * C0) * math.log2(2 * C0)
    rr_filter = C * rr_const * M(4 * C0) * math.log2(4 * C0)
    rr_units = (rr_linear + rr_filter) / M(flat_length)

    expected_candidates = 0.5 + cap * frobenius_filter_failure
    candidate_units = (1.0 + 2.0 * a) * C * M(precision) / M(flat_length)
    terminal_units = (
        tree_units
        + rr_units
        + expected_candidates * candidate_units
        + MISC_LONG_PRODUCTS
    )
    terminal_expected = (
        restart * terminal_units * long_extension_product_gates(flat_length, d)
    )

    # Re-lift all base-field branches only for surviving scalar candidates.
    coordinate_relifts = restart * expected_candidates * all_branches

    # Form lambda.x, the two Frobenius sketches, and filter leaves.  Full
    # extension multiplication is an overcharge for multiplication by F_127
    # scalars, but keeps the ledger simple.
    leaf_formation = (
        restart
        * C
        * (3.0 * a + 16.0)
        * precision
        * strict.sparse_extension_gate_envelope(d)
    )

    auxiliary = restart * (2.0**old.auxiliary_stress)
    preprocessing_log = strict.preprocessing_log2(v_expanded, m, r, restart)
    preprocessing = 2.0**preprocessing_log

    costs = {
        "initial_base_field_branch_lifts": math.log2(branch_expected),
        "valuation_safe_filter_series": math.log2(rational_filter_series),
        "terminal_scalar_recombination": math.log2(terminal_expected),
        "on_demand_coordinate_relifts": math.log2(coordinate_relifts),
        "extension_leaf_formation": math.log2(leaf_formation),
        "auxiliary_factorization_and_verification": math.log2(auxiliary),
        "preprocessing": preprocessing_log,
    }
    total = log2sum(*costs.values())
    costs.update(
        {
            "total": total,
            "target": float(target),
            "margin": float(target) - total,
            "strict_bundle_baseline": old.total_stress,
            "saving_vs_strict_bundle": old.total_stress - total,
        }
    )

    return Row(
        level=level,
        target=target,
        a=a,
        r=r,
        extension_degree=d,
        C=C,
        C0=C0,
        precision=precision,
        root_cap=cap,
        probabilities={
            "eligible_nonsingular_projective_root_lower": p_root,
            "root_and_cap_lower": p_root_cap,
            "chart": p_chart,
            "seed_disjoint": p_seed,
            "separator_failure_upper": sep_failure,
            "valuation_denominator_failure_upper": denominator_failure,
            "frobenius_descent_filter_failure_upper": frobenius_filter_failure,
            "one_trial_success_lower": math.exp(log_success),
            "restart_factor": restart,
        },
        multiplication={
            "relaxed_product_R_over_M": relaxed / M(precision),
            "long_extension_product_log2_gates": math.log2(
                long_extension_product_gates(flat_length, d)
            ),
        },
        terminal={
            "value_families": value_families,
            "tree_long_product_units": tree_units,
            "rational_reconstruction_long_product_units": rr_units,
            "candidate_long_product_units_each": candidate_units,
            "expected_candidates_per_trial_upper": expected_candidates,
            "misc_long_product_units": MISC_LONG_PRODUCTS,
            "total_long_product_units_per_trial": terminal_units,
        },
        cost_log2_gates=costs,
        sensitivity={},
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("QRUOV_split_lift_audit.json"),
    )
    args = parser.parse_args()

    rows = [attack_row(*params) for params in PARAMS]
    for row, params in zip(rows, PARAMS):
        for rr in (10.0, 16.0, 32.0, 64.0, 128.0):
            stressed = attack_row(*params, rr_const=rr)
            row.sensitivity[f"total_with_RR_constant_{int(rr)}"] = (
                stressed.cost_log2_gates["total"]
            )

    payload = {
        "description": (
            "Split-before-lift QR-UOV direct-forgery stress ledger with "
            "Frobenius descent and a valuation-safe rational oil filter."
        ),
        "status": (
            "Candidate finite gate audit under the strict bundle's asymptotic "
            "multiplication convention; not a measured circuit."
        ),
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print("QR-UOV split-before-lift audit")
    print("=" * 72)
    for row in rows:
        c = row.cost_log2_gates
        print(
            f"Level {row.level}: a={row.a}, d={row.extension_degree}, "
            f"cost={c['total']:.6f}, target={row.target}, "
            f"margin={c['margin']:+.6f}, saving={c['saving_vs_strict_bundle']:.6f}"
        )
        print(
            f"  branch={c['initial_base_field_branch_lifts']:.3f}, "
            f"terminal={c['terminal_scalar_recombination']:.3f}, "
            f"coordinates={c['on_demand_coordinate_relifts']:.3f}, "
            f"restart={row.probabilities['restart_factor']:.6f}"
        )
        print(
            "  RR sensitivity: "
            + ", ".join(
                f"{name.rsplit('_', 1)[-1]}->{value:.3f}"
                for name, value in row.sensitivity.items()
            )
        )


if __name__ == "__main__":
    main()
