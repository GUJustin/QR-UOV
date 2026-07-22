#!/usr/bin/env python3
"""Reproduce every new numerical bound in the corrected QR-UOV paper.

The script covers the random-projective-slice theorem and every numerical
cost column retained in the corrected manuscript.  Removed sensitivity
columns are intentionally not reproduced.

All logarithms are base two.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, getcontext
import math

getcontext().prec = 120

Q_INT = 127**3
Q = Decimal(Q_INT)
LOG2_Q = math.log2(Q_INT)
EPSILON = 0.2


@dataclass(frozen=True)
class Params:
    level: str
    V: int
    O: int
    m: int
    target: int
    solver_degree: int  # [L:K]


PARAMS = (
    Params("I", 52, 18, 54, 143, 10),
    Params("III", 76, 26, 78, 207, 13),
    Params("V", 102, 35, 105, 272, 16),
)


def log2_decimal(x: Decimal) -> float:
    if x <= 0:
        raise ValueError("logarithm argument must be positive")
    return float(x.ln() / Decimal(2).ln())


def rank_deficiency_log2(rows: int, cols: int) -> float:
    """log2 Pr[uniform rows-by-cols matrix has rank < cols]."""
    if rows < cols:
        return 0.0
    full = Decimal(1)
    for j in range(cols):
        full *= Decimal(1) - Q ** Decimal(j - rows)
    return log2_decimal(Decimal(1) - full)


def dependence_log2_upper(form_dim: int, m: int) -> float:
    """Log of the union-bound upper bound for m dependent uniform vectors."""
    numerator = Q**m - Decimal(1)
    bound = Q ** Decimal(-form_dim) * numerator / (Q - Decimal(1))
    return log2_decimal(bound)


def log2sum(log_values: list[float]) -> float:
    top = max(log_values)
    return top + math.log2(sum(2.0 ** (x - top) for x in log_values))


def fast_direct_indicator(V: int, epsilon: float = EPSILON) -> float:
    """Specialization of van der Hoeven--Lecerf Corollary 5.5."""
    return (
        V
        + (1.0 + epsilon) * (V - 1)
        + math.log2(3.0 * math.log2(127.0))
    )


def arbitrary_core_indicators(
    variables: int, base_degree: int, extension_degree: int, omega: float
) -> tuple[float, float]:
    """Direct/uniform indicators from Theorem (arbitrary-core lifting)."""
    coefficient_log = (
        math.log2(variables**omega + variables**2)
        + (2 * variables - 1)
        + math.log2(variables + 2)
    )
    field_bitlength = extension_degree * base_degree * math.log2(127.0)
    direct = coefficient_log + math.log2(field_bitlength)
    return direct, direct + math.log2(direct)


def row(p: Params) -> dict[str, float | int | str]:
    V, O, m = p.V, p.O, p.m
    e = m - V

    # Off-oil first-jet incidence.
    c_det = (e + 1) * O
    c_off = e * (O + 1) + 1
    log2_cumulative_degree = (
        math.log2(V)
        + m * math.log2(3.0)
        + c_det * math.log2(2.0 * V)
    )
    log2_eta_off = log2_cumulative_degree - c_off * LOG2_Q

    # No regular geometric oil direction implies all standard directions bad.
    log2_eta_oil = O * rank_deficiency_log2(m, V)

    # Linear dependence among the m central forms.
    form_dim = V * (V + 1) // 2 + V * O
    log2_eta_dep = dependence_log2_upper(form_dim, m)

    log2_key_bad = log2sum([log2_eta_off, log2_eta_oil, log2_eta_dep])

    # Full discriminant after quadratic slice substitution and linear mixing.
    partial_degree = (V + 1) * (2 ** (V - 1))
    joint_degree = 3 * V * partial_degree
    log2_joint_degree = math.log2(joint_degree)
    log2_field_size = p.solver_degree * LOG2_Q
    log2_delta = log2_joint_degree - log2_field_size

    direct = fast_direct_indicator(V)
    uniform = direct + math.log2(direct)

    fallback_extension = math.ceil((3.0 + 2.0 * V) / LOG2_Q)
    fallback_direct, fallback_uniform = arbitrary_core_indicators(
        V, 3, fallback_extension, math.log2(7.0)
    )
    _, fallback_cubic_uniform = arbitrary_core_indicators(
        V, 3, fallback_extension, 3.0
    )

    residual = m - 3
    residual_extension = math.ceil(
        (3.0 + 2.0 * residual) / math.log2(127.0)
    )
    pseudo_direct, pseudo_uniform = arbitrary_core_indicators(
        residual, 1, residual_extension, math.log2(7.0)
    )

    return {
        "level": p.level,
        "V": V,
        "O": O,
        "m": m,
        "e": e,
        "c_det": c_det,
        "c_off": c_off,
        "form_dim": form_dim,
        "log2_cumulative_degree": log2_cumulative_degree,
        "log2_eta_off": log2_eta_off,
        "log2_eta_oil": log2_eta_oil,
        "log2_eta_dep": log2_eta_dep,
        "log2_key_bad": log2_key_bad,
        "joint_degree": joint_degree,
        "log2_joint_degree": log2_joint_degree,
        "log2_field_size": log2_field_size,
        "log2_delta": log2_delta,
        "log2_delta_two": 2.0 * log2_delta,
        "fast_direct": direct,
        "fast_uniform": uniform,
        "uniform_margin": p.target - uniform,
        "fallback_extension": fallback_extension,
        "fallback_direct": fallback_direct,
        "fallback_uniform": fallback_uniform,
        "fallback_cubic_uniform": fallback_cubic_uniform,
        "residual": residual,
        "residual_extension": residual_extension,
        "pseudo_direct": pseudo_direct,
        "pseudo_uniform": pseudo_uniform,
        "target": p.target,
    }


def main() -> None:
    rows = [row(p) for p in PARAMS]
    print(f"Q = {Q_INT}")
    print(f"log2(Q) = {LOG2_Q:.12f}")
    print(f"fast-solver epsilon = {EPSILON}\n")

    print("KEY-LEVEL BAD-SET BOUNDS")
    print(
        "level  c_off  log2(degree factor)  log2 eta_off  "
        "log2 eta_oil  log2 eta_dep  log2 total"
    )
    for r in rows:
        print(
            f"{r['level']:>5}  {r['c_off']:>5}  "
            f"{r['log2_cumulative_degree']:>19.3f}  "
            f"{r['log2_eta_off']:>12.3f}  "
            f"{r['log2_eta_oil']:>12.3f}  "
            f"{r['log2_eta_dep']:>12.3f}  "
            f"{r['log2_key_bad']:>10.3f}"
        )

    print("\nJOINT SLICE/OUTPUT REGULARIZATION")
    print("level  [L:K]  log2 degree  log2 |L|  log2 delta  two trials")
    for p, r in zip(PARAMS, rows):
        print(
            f"{r['level']:>5}  {p.solver_degree:>5}  "
            f"{r['log2_joint_degree']:>11.3f}  "
            f"{r['log2_field_size']:>8.3f}  "
            f"{r['log2_delta']:>10.3f}  "
            f"{r['log2_delta_two']:>10.3f}"
        )

    print("\nFAST COST LEDGER")
    print("level  direct  uniform TlogT  uniform margin  target")
    for r in rows:
        print(
            f"{r['level']:>5}  {r['fast_direct']:>6.2f}  "
            f"{r['fast_uniform']:>13.2f}  "
            f"{r['uniform_margin']:>14.2f}  {r['target']:>6}"
        )

    print("\nARBITRARY-CORE FALLBACK LEDGER")
    print("level  [K':K]  Strassen direct/uniform  cubic uniform  target")
    for r in rows:
        print(
            f"{r['level']:>5}  {r['fallback_extension']:>6}  "
            f"{r['fallback_direct']:>8.2f}/{r['fallback_uniform']:<8.2f}  "
            f"{r['fallback_cubic_uniform']:>13.2f}  {r['target']:>6}"
        )

    print("\nPSEUDO-OIL TERMINAL-SOLVE LEDGER")
    print("level  residual  [K':F_127]  Strassen direct/uniform  target")
    for r in rows:
        print(
            f"{r['level']:>5}  {r['residual']:>8}  "
            f"{r['residual_extension']:>10}  "
            f"{r['pseudo_direct']:>8.2f}/{r['pseudo_uniform']:<8.2f}  "
            f"{r['target']:>6}"
        )


if __name__ == "__main__":
    main()
