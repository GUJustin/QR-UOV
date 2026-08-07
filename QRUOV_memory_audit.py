#!/usr/bin/env python3
"""Reconstruct time and memory scales for QR-UOV attack estimates.

The cited papers mostly report arithmetic or gate counts.  This script derives
the dimensions of their dominant linear-algebra objects from the published
formulas and reports three storage models:

* vector: three packed field vectors, a conservative Wiedemann working state;
* sparse: one stored sparse matrix with explicit column indices;
* dense: one dense square matrix.

For PXL, the script reports the packed size of the dense alpha-by-alpha state
appearing in the published alpha^omega elimination term.  For Just Guess, it
reports the packed public quadratic coefficients plus two change-of-basis
matrices; its depth-first search stack is lower order.

These are order-of-magnitude reconstructions, not measured peak resident sets.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass


Q = 127
LOGQ = math.log2(Q)
FIELD_BITS = 7
OMEGA = 2.37


@dataclass
class Row:
    level: str
    method: str
    time_log2: float
    object_log2: float | None
    vector_log2_bytes: float | None
    sparse_log2_bytes: float | None
    dense_log2_bytes: float | None
    details: dict[str, int | float | str]


def semiregular_degree(n: int, m: int, guessed: int) -> int:
    """First non-positive coefficient of (1-z^2)^m/(1-z)^(n-k+1)."""
    denominator = n - guessed + 1
    for degree in range(1, 1000):
        coefficient = 0
        for j in range(degree // 2 + 1):
            if j > m:
                break
            coefficient += (
                (-1) ** j
                * math.comb(m, j)
                * math.comb(denominator + degree - 2 * j - 1, degree - 2 * j)
            )
        if coefficient <= 0:
            return degree
    raise RuntimeError("degree search did not terminate")


def field_gate_log2(field_log2: float) -> float:
    return math.log2(2 * field_log2 * field_log2 + field_log2)


def linear_algebra_memory(
    log_n: float, rho: int, field_bits: int
) -> tuple[float, float, float]:
    """Three vectors, stored sparse matrix, and dense square matrix in bytes."""
    field_bytes = math.ceil(field_bits / 8)
    index_bytes = math.ceil(log_n / 8)
    vector = log_n + math.log2(3 * field_bits / 8)
    sparse = log_n + math.log2(rho) + math.log2(field_bytes + index_bytes)
    dense = 2 * log_n + math.log2(field_bytes)
    return vector, sparse, dense


def best_wxl(n: int, m: int, field_log2: float = LOGQ) -> dict[str, float | int]:
    candidates = []
    # For a square affine system, the k=0 semi-regular series has no
    # non-positive coefficient under this convention.  Hybrid WXL therefore
    # starts with one guessed variable in that case.
    first_guess = max(0, n - m)
    if n == m:
        first_guess = 1
    for guessed in range(first_guess, n + 1):
        degree = semiregular_degree(n, m, guessed)
        variables = n - guessed
        dimension = math.comb(degree + variables, degree)
        rho = math.comb(variables + 2, 2)
        time = (
            guessed * field_log2
            + math.log2(3)
            + math.log2(rho)
            + 2 * math.log2(dimension)
            + field_gate_log2(field_log2)
        )
        candidates.append((time, guessed, degree, dimension, rho, variables))
    time, guessed, degree, dimension, rho, variables = min(candidates)
    return {
        "time": time,
        "guessed": guessed,
        "degree": degree,
        "dimension": dimension,
        "rho": rho,
        "variables_after_guessing": variables,
    }


def pxl_alpha(n: int, m: int, guessed: int, degree: int) -> int:
    exponent = m - (n - guessed)
    total = 0
    for d in range(degree + 1):
        coefficient = 0
        for j in range(max(0, d - m), min(exponent, d) + 1):
            coefficient += (-1) ** j * math.comb(exponent, j) * math.comb(m, d - j)
        total += max(coefficient, 0)
    return total


def best_pxl(n: int, m: int) -> dict[str, float | int]:
    candidates = []
    for guessed in range(1, n + 1):
        degree = semiregular_degree(n, m, guessed)
        alpha = pxl_alpha(n, m, guessed, degree)
        left = math.comb(n - guessed + degree, degree)
        right = math.comb(n + degree, degree)
        guessed_monomials = math.comb(guessed + degree, degree)
        terms = [
            2 * math.log2(guessed)
            + math.log2(alpha)
            + math.log2(left)
            + math.log2(right),
            guessed * LOGQ
            + 2 * math.log2(alpha)
            + math.log2(guessed_monomials),
            OMEGA * math.log2(alpha),
        ]
        top = max(terms)
        arithmetic = top + math.log2(sum(2 ** (term - top) for term in terms))
        time = arithmetic + field_gate_log2(LOGQ)
        candidates.append((time, guessed, degree, alpha, left, right, guessed_monomials))
    time, guessed, degree, alpha, left, right, guessed_monomials = min(candidates)
    return {
        "time": time,
        "guessed": guessed,
        "degree": degree,
        "alpha": alpha,
        "left_monomials": left,
        "right_monomials": right,
        "guessed_monomials": guessed_monomials,
    }


def pxl_row(level: str, method: str, n: int) -> Row:
    estimate = best_pxl(n, n)
    log_alpha = math.log2(int(estimate["alpha"]))
    dense_state = 2 * log_alpha + math.log2(FIELD_BITS / 8)
    return Row(
        level,
        method,
        float(estimate["time"]),
        log_alpha,
        None,
        None,
        dense_state,
        {**estimate, "variables": n, "memory_object": "packed alpha-by-alpha state"},
    )


def wxl_row(level: str, method: str, n: int, m: int, time_override: float | None = None) -> Row:
    estimate = best_wxl(n, m)
    log_n = math.log2(int(estimate["dimension"]))
    vector, sparse, dense = linear_algebra_memory(log_n, int(estimate["rho"]), FIELD_BITS)
    return Row(
        level,
        method,
        float(time_override if time_override is not None else estimate["time"]),
        log_n,
        vector,
        sparse,
        dense,
        {**estimate, "equations": m, "variables": n},
    )


def fixed_wxl_row(level: str, n: int, guessed: int) -> Row:
    degree = semiregular_degree(n, n, guessed)
    variables = n - guessed
    dimension = math.comb(degree + variables, degree)
    rho = math.comb(variables + 2, 2)
    time = (
        guessed * LOGQ
        + math.log2(3)
        + math.log2(rho)
        + 2 * math.log2(dimension)
        + field_gate_log2(LOGQ)
    )
    # Count three Wiedemann vectors and the reduced quadratic system using one
    # 64-bit machine word per F_127 element.  Matrix entries are regenerated.
    polynomial_coefficients = n * math.comb(variables + 2, 2)
    working_bytes = 8 * (3 * dimension + polynomial_coefficients)
    return Row(
        level,
        "This work, low-memory core WXL",
        time,
        math.log2(dimension),
        math.log2(working_bytes),
        None,
        None,
        {
            "core_variables": n,
            "guessed": guessed,
            "variables_after_guessing": variables,
            "degree": degree,
            "dimension": dimension,
            "rho": rho,
            "stored_field_elements": 3 * dimension + polynomial_coefficients,
            "memory_model": "three vectors and reduced quadratics as 64-bit words; matrix regenerated",
            "time_caveat": "matrix-generation overhead is not included in the published WXL field-gate formula",
        },
    )


def reconciliation_row(level: str, n: int, m: int, published_time: float) -> Row:
    field_log2 = 3 * LOGQ
    estimate = best_wxl(n, m, field_log2)
    log_n = math.log2(int(estimate["dimension"]))
    vector, sparse, dense = linear_algebra_memory(log_n, int(estimate["rho"]), 3 * FIELD_BITS)
    return Row(
        level,
        "Reconciliation WXL",
        published_time,
        log_n,
        vector,
        sparse,
        dense,
        {**estimate, "equations": m, "variables": n, "field": "F_(127^3)"},
    )


def rectangular_minrank_row(
    level: str,
    n: int,
    m: int,
    rank: int,
    dimension: int,
    degree: int,
    published_time: float,
) -> Row:
    variables = n - dimension + 1
    macaulay = math.comb(variables + degree - 1, degree) * math.comb(m, rank)
    rho = (rank + 1) * variables
    log_n = math.log2(macaulay)
    vector, sparse, dense = linear_algebra_memory(log_n, rho, 3 * FIELD_BITS)
    return Row(
        level,
        "Rectangular MinRank",
        published_time,
        log_n,
        vector,
        sparse,
        dense,
        {
            "uov_variables": n,
            "equations": m,
            "target_rank": rank,
            "solution_dimension": dimension,
            "macaulay_degree": degree,
            "remaining_variables": variables,
            "dimension": macaulay,
            "rho": rho,
            "field": "F_(127^3)",
        },
    )


def wedge_row(level: str, v: int, reduced_oil: int, published_time: float) -> Row:
    dimension = (
        math.comb(v + reduced_oil, v)
        * math.comb(v + reduced_oil + 1, v)
        // (v + 1)
    )
    rho = (v + 1) ** 2
    log_n = math.log2(dimension)
    vector, sparse, dense = linear_algebra_memory(log_n, rho, 3 * FIELD_BITS)
    return Row(
        level,
        "Wedge",
        published_time,
        log_n,
        vector,
        sparse,
        dense,
        {
            "vinegar_variables": v,
            "reduced_oil": reduced_oil,
            "dimension": dimension,
            "rho": rho,
            "field": "F_(127^3)",
            "time_unit": "field operations",
        },
    )


def just_guess_row(level: str, n: int, m: int, published_time: float) -> Row:
    coefficients = m * n * (n + 1) // 2 + n * n + m * m
    packed = math.log2(coefficients * FIELD_BITS / 8)
    machine_words = math.log2(coefficients * 8)
    return Row(
        level,
        "Just Guess",
        published_time,
        math.log2(coefficients),
        packed,
        None,
        machine_words,
        {
            "variables": n,
            "equations": m,
            "stored_field_elements": coefficients,
            "memory_object": "public quadratics plus S and T; DFS stack is lower order",
            "dense_log2_bytes_uses": "one 64-bit machine word per field element",
        },
    )


def rows() -> list[Row]:
    levels = {
        "I": {
            "full_n": 210,
            "m": 54,
            "core": 50,
            "pxl_n": 52,
            "hash": (49, 51, 158),
            "recon": (52, 54, 164),
            "rm": (70, 54, 52, 18, 20, 158),
            "wedge": (52, 13, 182),
            "jg": 258,
            "low_memory_guess": 27,
        },
        "III": {
            "full_n": 306,
            "m": 78,
            "core": 74,
            "pxl_n": 76,
            "hash": (72, 75, 217),
            "recon": (76, 78, 231),
            "rm": (102, 78, 76, 26, 29, 219),
            "wedge": (76, 17, 249),
            "jg": 399,
            "low_memory_guess": 45,
        },
        "V": {
            "full_n": 411,
            "m": 105,
            "core": 101,
            "pxl_n": 103,
            "hash": (99, 102, 285),
            "recon": (102, 105, 291),
            "rm": (137, 105, 102, 35, 35, 277),
            "wedge": (102, 22, 328),
            "jg": 561,
            "low_memory_guess": 70,
        },
    }
    output = []
    for level, p in levels.items():
        output.append(pxl_row(level, "Specification PXL", int(p["pxl_n"])))
        output.append(pxl_row(level, "This work, core plus PXL", int(p["core"])))
        h_n, h_m, h_time = p["hash"]
        output.append(wxl_row(level, "Updated Hashimoto WXL", h_n, h_m, h_time))
        output.append(wxl_row(level, "This work, core plus WXL", int(p["core"]), int(p["core"])))
        output.append(fixed_wxl_row(level, int(p["core"]), int(p["low_memory_guess"])))
        r_n, r_m, r_time = p["recon"]
        output.append(reconciliation_row(level, r_n, r_m, r_time))
        output.append(rectangular_minrank_row(level, *p["rm"]))
        output.append(wedge_row(level, *p["wedge"]))
        output.append(just_guess_row(level, int(p["full_n"]), int(p["m"]), float(p["jg"])))
    return output


def validate(result: list[Row]) -> None:
    expected_pxl = {"I": 150, "III": 211, "V": 279}
    expected_reconciliation = {"I": 164, "III": 231, "V": 291}
    expected_rm = {"I": 158, "III": 219, "V": 277}
    expected_wedge = {"I": 182, "III": 249, "V": 328}
    for row in result:
        if row.method == "Specification PXL":
            assert abs(row.time_log2 - expected_pxl[row.level]) < 0.6
        elif row.method == "Reconciliation WXL":
            assert row.time_log2 == expected_reconciliation[row.level]
        elif row.method == "Rectangular MinRank":
            assert row.time_log2 == expected_rm[row.level]
        elif row.method == "Wedge":
            assert row.time_log2 == expected_wedge[row.level]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=str)
    args = parser.parse_args()
    result = rows()
    validate(result)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump([asdict(row) for row in result], handle, indent=2)
            handle.write("\n")
    for row in result:
        fields = [
            f"{row.level:>3}",
            f"{row.method:<34}",
            f"time=2^{row.time_log2:.3f}",
        ]
        if row.object_log2 is not None:
            fields.append(f"log2(object)={row.object_log2:.3f}")
        if row.vector_log2_bytes is not None:
            fields.append(f"vector=2^{row.vector_log2_bytes:.3f} B")
        if row.sparse_log2_bytes is not None:
            fields.append(f"sparse=2^{row.sparse_log2_bytes:.3f} B")
        if row.dense_log2_bytes is not None:
            fields.append(f"dense=2^{row.dense_log2_bytes:.3f} B")
        print(" | ".join(fields))


if __name__ == "__main__":
    main()
