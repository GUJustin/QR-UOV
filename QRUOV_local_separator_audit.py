#!/usr/bin/env python3
"""Local-separator split-before-lift stress ledger.

This variant replaces global pairwise separation of all C branches by two
one-sided requirements relative to a canonical eligible target root:
  (i) lambda has the correct leading valuation on every branch at infinity;
  (ii) lambda separates the canonical target root from every other finite root.
It accepts only simple target factors. A single Frobenius sketch using lambda
itself then gives exact base-field descent on such factors.

The v2 correctness pass adds an explicit one-time charge for forming separator
leaves during coordinate recovery and records conservative root-probability and
online-multiplication sensitivities.
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
BASE_PATH = HERE / "QRUOV_split_lift_audit.py"
spec = importlib.util.spec_from_file_location("splitbase", BASE_PATH)
base = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = base
spec.loader.exec_module(base)
strict = base.strict
Q = base.Q
LOGQ = base.LOGQ
GMUL = base.GMUL
GADD = base.GADD

# The parent manuscript's already-proved projective-root lower bound. The
# sharper second-Bonferroni value remains the primary ledger value, but this
# constant gives an independent conservative fallback.
PARENT_ROOT_LOWER = 0.326336998767509

# Audited sparse irreducibles over F_127 found by exhaustive search among
# monic polynomials with at most three nonleading +/-1 terms.
EXTRA_MODULI = {
    7: {5: 1, 1: 1, 0: 1},
    8: {4: 1, 0: -1},
    9: {2: 1, 0: 1},
    10: {6: -1, 4: 1, 0: 1},
    11: {5: -1, 0: 1},
    12: {3: -1, 1: 1, 0: 1},
    13: {4: 1, 1: -1, 0: 1},
    14: {11: 1, 1: 1, 0: -1},
}
strict.base.SPARSE_MODULI.update(EXTRA_MODULI)


@dataclass
class LocalRow:
    level: str
    target: int
    a: int
    r: int
    extension_degree: int
    modulus: str
    C: int
    C0: int
    precision: int
    probabilities: dict[str, float]
    terminal: dict[str, float]
    cost_log2_gates: dict[str, float]
    sensitivity: dict[str, object]


def attack_row(
    params,
    d: int,
    *,
    rr_const: float = base.RR_CONST,
    root_lower: float | None = None,
    relaxed_multiplier: float = 1.0,
) -> LocalRow:
    level, v_expanded, m, V, O, target = params
    a = m - 4
    r = m - 3
    C = 1 << a
    C0 = (1 << (a - 1)) * (a + 2)
    precision = 4 * C0 + 1
    flat_length = 4 * C * precision
    K = float(Q**d)
    if d not in strict.base.SPARSE_MODULI:
        raise KeyError(d)
    if 2.0 * C >= K:
        raise ValueError("local separator field too small")
    if relaxed_multiplier <= 0:
        raise ValueError("relaxed_multiplier must be positive")

    old = strict.direct_row(level, v_expanded, m, V, O, target)
    p_root = base.projective_root_lower(r) if root_lower is None else root_lower
    p_chart = 1.0 - 1.0 / Q
    p_seed = strict.seed_disjoint_probability_lower(v_expanded, m)

    # If I branches escape and F remain finite, I+F=C. There are at most I
    # valuation hyperplanes and F-1 collisions against P_star.
    lambda_failure = (C - 1.0) / K

    # Sample c uniformly from K' minus the known forbidden values
    # {-lambda(xi_s)} at the start fiber. Conditional on this choice, at most
    # F<=C target values remain forbidden.
    denominator_failure = C / (K - C)
    log_success = (
        math.log(p_seed)
        + math.log(p_root)
        + math.log(p_chart)
        + base.log_success_from_failure(lambda_failure)
        + base.log_success_from_failure(denominator_failure)
    )
    restart = math.exp(-log_success)

    h = a * (a + 1) // 2
    relaxed = relaxed_multiplier * base.relaxed_upper(precision)
    one_branch = h * relaxed * GMUL + (a * h + a * a + 8 * a) * precision * (GMUL + GADD)
    all_branches = C * one_branch
    branch_expected = restart * all_branches

    rational_filter_series = (
        restart
        * C
        * base.FILTER_INVERSE_CONST
        * base.M(precision)
        * strict.sparse_extension_gate_envelope(d)
    )

    # Families: b1=lambda^(q^{-1}).x and psi=H0/(c+lambda.x).
    value_families = 2
    tree_units = (1.0 + 2.0 * value_families) * a / 2.0

    # Rational reconstructions: q and b1 have degree <=C0; psi graph <=2C0.
    rr_linear = 2.0 * C * rr_const * base.M(2 * C0) * math.log2(2 * C0)
    rr_filter = C * rr_const * base.M(4 * C0) * math.log2(4 * C0)
    rr_units = (rr_linear + rr_filter) / base.M(flat_length)

    # One coordinate product-tree evaluation is performed only on the final
    # successful trial. Its long-product work is not multiplied by restart.
    candidate_units = (1.0 + 2.0 * a) * C * base.M(precision) / base.M(flat_length)
    fixed_terminal_units = tree_units + rr_units + base.MISC_LONG_PRODUCTS
    long_gates = base.long_extension_product_gates(flat_length, d)
    terminal_expected = restart * fixed_terminal_units * long_gates + candidate_units * long_gates
    terminal_units = fixed_terminal_units + candidate_units / restart

    # Exactly one complete base-field branch re-lift is needed over the whole
    # Las Vegas execution.
    coordinate_relifts = all_branches

    ext_scalar_gate = strict.sparse_extension_gate_envelope(d)
    scalar_leaf_formation = restart * C * (2.0 * a + 16.0) * precision * ext_scalar_gate

    # v1 omitted this small cost: the coordinate rerun must form u_s=lambda.x_s
    # again before its product tree. Coordinate-series leaves themselves are
    # base-field embeddings and need no extension multiplication.
    coordinate_leaf_formation = C * a * precision * ext_scalar_gate

    # Retain the parent auxiliary estimate unchanged: it uses the larger old
    # field and is therefore a conservative overcharge for factorization and
    # exact public verification in this variant.
    auxiliary = restart * (2.0**old.auxiliary_stress)
    preprocessing_log = strict.preprocessing_log2(v_expanded, m, r, restart)
    costs = {
        "initial_base_field_branch_lifts": math.log2(branch_expected),
        "valuation_safe_filter_series": math.log2(rational_filter_series),
        "terminal_scalar_recombination": math.log2(terminal_expected),
        "on_demand_coordinate_relifts": math.log2(coordinate_relifts),
        "extension_leaf_formation": math.log2(scalar_leaf_formation),
        "coordinate_relift_leaf_formation": math.log2(coordinate_leaf_formation),
        "auxiliary_factorization_and_verification": math.log2(auxiliary),
        "preprocessing": preprocessing_log,
    }
    total = base.log2sum(*costs.values())
    costs.update(
        total=total,
        target=float(target),
        margin=float(target) - total,
        split_global_baseline=base.attack_row(*params).cost_log2_gates["total"],
        saving_vs_global_split=base.attack_row(*params).cost_log2_gates["total"] - total,
        strict_bundle_baseline=old.total_stress,
        saving_vs_strict_bundle=old.total_stress - total,
    )
    return LocalRow(
        level,
        target,
        a,
        r,
        d,
        strict.base.modulus_string(d),
        C,
        C0,
        precision,
        {
            "eligible_nonsingular_projective_root_lower": p_root,
            "chart": p_chart,
            "seed_disjoint": p_seed,
            "valuation_or_target_collision_failure_upper": lambda_failure,
            "valuation_denominator_failure_upper": denominator_failure,
            "one_trial_success_lower": math.exp(log_success),
            "restart_factor": restart,
        },
        {
            "value_families": value_families,
            "tree_long_product_units": tree_units,
            "rational_reconstruction_long_product_units": rr_units,
            "candidate_long_product_units_each": candidate_units,
            "coordinate_candidates_per_success": 1.0,
            "misc_long_product_units": base.MISC_LONG_PRODUCTS,
            "total_long_product_units_per_trial": terminal_units,
            "relaxed_product_multiplier": relaxed_multiplier,
        },
        costs,
        {},
    )


def candidate_rows(params, *, rr_const=base.RR_CONST, root_lower=None, relaxed_multiplier=1.0):
    rows = []
    oldd = strict.direct_row(*params).working_degree
    for d in range(2, oldd + 1):
        if d not in strict.base.SPARSE_MODULI:
            continue
        try:
            rows.append(
                attack_row(
                    params,
                    d,
                    rr_const=rr_const,
                    root_lower=root_lower,
                    relaxed_multiplier=relaxed_multiplier,
                )
            )
        except ValueError:
            pass
    return rows


def best_row(params, **kwargs):
    return min(candidate_rows(params, **kwargs), key=lambda x: x.cost_log2_gates["total"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", type=Path, default=Path("QRUOV_local_separator_audit.json"))
    args = ap.parse_args()

    rows = []
    for params in base.PARAMS:
        candidates = candidate_rows(params)
        row = min(candidates, key=lambda x: x.cost_log2_gates["total"])

        for rr in (10.0, 16.0, 32.0, 64.0, 128.0, 256.0):
            stressed = best_row(params, rr_const=rr)
            row.sensitivity[f"total_with_RR_constant_{int(rr)}"] = {
                "extension_degree": stressed.extension_degree,
                "total": stressed.cost_log2_gates["total"],
            }

        row.sensitivity["degree_candidates"] = {
            str(x.extension_degree): x.cost_log2_gates["total"] for x in candidates
        }

        conservative = best_row(params, root_lower=PARENT_ROOT_LOWER)
        row.sensitivity["parent_root_probability_fallback"] = {
            "root_lower": PARENT_ROOT_LOWER,
            "extension_degree": conservative.extension_degree,
            "total": conservative.cost_log2_gates["total"],
            "margin": conservative.cost_log2_gates["margin"],
            "restart_factor": conservative.probabilities["restart_factor"],
        }

        online = {}
        for mult in (2.0, 4.0, 8.0, 16.0):
            stressed = best_row(params, relaxed_multiplier=mult)
            online[str(int(mult))] = {
                "extension_degree": stressed.extension_degree,
                "total": stressed.cost_log2_gates["total"],
                "margin": stressed.cost_log2_gates["margin"],
            }
        row.sensitivity["online_product_envelope_multiplier"] = online
        rows.append(row)

    certs = []
    for d in sorted({r.extension_degree for r in rows} | {8, 9, 11, 12, 15, 16}):
        certs.append(strict.base.rabin_irreducibility_certificate(d))

    payload = {
        "description": "Local-separator split-before-lift QR-UOV stress ledger with exact simple-root Frobenius descent.",
        "version": "v2 correctness pass",
        "status": "Candidate finite gate audit; proof still depends on exact applicability of the characteristic-polynomial degree theorem, the rational-filter graph bound, and a fully explicit online-multiplication schedule.",
        "rows": [asdict(r) for r in rows],
        "sparse_modulus_certificates": certs,
    }
    args.json.write_text(json.dumps(payload, indent=2) + "\n")

    for r in rows:
        c = r.cost_log2_gates
        fallback = r.sensitivity["parent_root_probability_fallback"]
        print(
            f"Level {r.level}: d={r.extension_degree}, cost={c['total']:.9f}, "
            f"margin={c['margin']:+.9f}, global-saving={c['saving_vs_global_split']:.6f}"
        )
        print(
            f"  terminal={c['terminal_scalar_recombination']:.3f}, "
            f"restart={r.probabilities['restart_factor']:.6f}, modulus={r.modulus}"
        )
        print(
            f"  parent-root fallback: d={fallback['extension_degree']}, "
            f"total={fallback['total']:.9f}, margin={fallback['margin']:+.9f}"
        )
        print("  degree candidates:", r.sensitivity["degree_candidates"])
        print("  online multipliers:", r.sensitivity["online_product_envelope_multiplier"])


if __name__ == "__main__":
    main()
