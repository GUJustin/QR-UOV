#!/usr/bin/env python3
"""Independent arithmetic audit for the QR-UOV paper.

Standard-library only.  Recomputes the paper's probability, field-size,
cost-indicator, seeded-transfer, and memory tables from the parameter tuples.
It is not an implementation of the polynomial-system solvers.
"""
from __future__ import annotations

from decimal import Decimal, getcontext
from math import ceil, log2

getcontext().prec = 120
D2 = Decimal(2)
LN2 = D2.ln()
q = 127
ell = 3
Q = q ** ell
LOGQ = log2(Q)
EIB = 2 ** 60

PARAMS = [
    dict(level="I", V=52, m=54, O=18, target=143, published=150,
         lam=128, tau1=4267, tau2=2916, rdir=3, bsign=22),
    dict(level="III", V=76, m=78, O=26, target=207, published=211,
         lam=192, tau1=9020, tau2=6123, rdir=4, bsign=31),
    dict(level="V", V=102, m=105, O=35, target=272, published=277,
         lam=256, tau1=16144, tau2=11018, rdir=4, bsign=40),
]


def dec_log2(x: Decimal) -> float:
    if x <= 0:
        raise ValueError("log of nonpositive value")
    return float(x.ln() / LN2)


def log2sum(a: float, b: float) -> float:
    m = max(a, b)
    return m + log2(2 ** (a - m) + 2 ** (b - m))


def rank_failure(base: int, rows: int, cols: int) -> Decimal:
    """Probability that a uniform rows x cols matrix (rows>=cols) is rank deficient."""
    B = Decimal(base)
    full = Decimal(1)
    for j in range(cols):
        full *= Decimal(1) - B ** Decimal(j - rows)
    return Decimal(1) - full


def square_rank_probability(base: int, n: int, r: int) -> Decimal:
    B = Decimal(base)
    value = B ** Decimal(-n * n)
    for i in range(r):
        value *= (B ** n - B ** i) ** 2 / (B ** r - B ** i)
    return value


def gaussian_binomial(n: int, k: int, base: int) -> Decimal:
    B = Decimal(base)
    out = Decimal(1)
    for i in range(k):
        out *= (B ** (n - i) - 1) / (B ** (k - i) - 1)
    return out


def vdhl_log_bound(V: int, epsilon: float) -> float:
    eta = epsilon / 3.0
    first = max(
        1.0,
        log2((1 + eta) * (8 / eta + 1) * V) + eta / 8.0,
    ) + log2(2 ** V - 1)
    second = log2(100) + 3 * V
    return max(first, second)


def minimum_s(V: int, epsilon: float) -> int:
    # Strict inequality Q^s > B.  None of the audited values lies on equality.
    return int(ceil(vdhl_log_bound(V, epsilon) / LOGQ))


def alpha_beta(V: int, s: int) -> tuple[Decimal, Decimal]:
    D = Decimal(2) ** V
    field = Decimal(Q) ** s
    alpha = (Decimal(5) * D ** 2 + Decimal(26) * D ** 3) / field
    beta = Decimal(3 * V * (V + 1)) * D / Decimal(2) / field
    return alpha, beta


def fast_indicators(V: int, epsilon: float) -> tuple[float, float]:
    C = V + (1 + epsilon) * (V - 1) + log2(ell * log2(q))
    return C, C + log2(C)


def epsilon_threshold(V: int, target: int) -> float:
    lo, hi = 0.0, 2.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if fast_indicators(V, mid)[1] < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def pair_signing_bound(m: int) -> tuple[Decimal, Decimal, Decimal, Decimal]:
    pi = square_rank_probability(q, m, m)
    rho_m1 = square_rank_probability(q, m, m - 1)
    sigma = Decimal(1) - pi
    chi = Decimal(q) ** Decimal(-(m - 1)) + (
        Decimal(1) - Decimal(q) ** Decimal(-(m - 1))
    ) / Decimal(q)
    zeta = rho_m1 / Decimal(q) * chi + (Decimal(1) - pi - rho_m1) * (Decimal(1) - pi)
    return pi, sigma, rho_m1, zeta


def ideal_cipher_distance(T: int) -> Decimal:
    space = Decimal(2) ** 128
    no_collision = Decimal(1)
    for j in range(T):
        no_collision *= Decimal(1) - Decimal(j) / space
    return Decimal(1) - no_collision


def assert_close(name: str, got: float, want: float, tol: float = 5e-4) -> None:
    if abs(got - want) > tol:
        raise AssertionError(f"{name}: got {got}, expected {want} +/- {tol}")


print(f"Q = {Q}; log2(Q) = {LOGQ:.12f}\n")

print("GRAPH-DIRECTION AND SIGNING LISTS")
for p in PARAMS:
    V, m, O = p["V"], p["m"], p["O"]
    tau = rank_failure(Q, m, V)
    one_log = dec_log2(tau)
    degree = V * (m - V)
    per_direction = degree + 1
    compact_systems = p["rdir"] * per_direction
    pi, sigma, _, zeta = pair_signing_bound(m)
    basis_all_log = V * dec_log2(sigma)
    basis_short_log = p["bsign"] * dec_log2(sigma)
    pair_log = (V // 2) * dec_log2(zeta)
    compact_rank_log = p["rdir"] * one_log
    compact_union_log = dec_log2(
        (Decimal(2) ** Decimal(compact_rank_log))
        + (Decimal(2) ** Decimal(basis_short_log))
    )
    print(
        f"{p['level']:>3}: condenser degree/list={degree}/{per_direction}; "
        f"compact systems={compact_systems}; all-direction systems={O*per_direction}\n"
        f"     log2 rank fail: one={one_log:.6f}, compact={compact_rank_log:.6f}, all={O*one_log:.6f}\n"
        f"     pi_m={float(pi):.15f}, E[basis tests]<={float(1/pi):.12f}\n"
        f"     log2 sign fail: short={basis_short_log:.6f}, all basis={basis_all_log:.6f}, pair fallback={pair_log:.6f}\n"
        f"     log2 compact union={compact_union_log:.6f}"
    )

# Frozen headline assertions.
expected_rank_one = [-62.898161, -62.898161, -83.864216]
expected_rank_all = [-1132.166907, -1635.352198, -2935.247544]
expected_basis_short = [-153.502126, -216.298451, -279.094775]
expected_pair = [-544.820, -796.276, -1068.687]
for p, e1, ea, es, ep in zip(PARAMS, expected_rank_one, expected_rank_all, expected_basis_short, expected_pair):
    tau = rank_failure(Q, p["m"], p["V"])
    pi, sigma, _, zeta = pair_signing_bound(p["m"])
    assert_close(f"{p['level']} one-direction", dec_log2(tau), e1, 6e-4)
    assert_close(f"{p['level']} all-directions", p["O"]*dec_log2(tau), ea, 6e-4)
    assert_close(f"{p['level']} short signing", p["bsign"]*dec_log2(sigma), es, 6e-4)
    assert_close(f"{p['level']} pair signing", (p["V"]//2)*dec_log2(zeta), ep, 1e-3)

print("\nFAST SOLVER FIELD AND ONE-INVOCATION ACCOUNTING (epsilon=0.3)")
expected_min_s = [8, 12, 15]
expected_min_success = [0.992337625602, 0.999997944684, 0.927725586161]
expected_safe_success = [0.999999996259, 0.999999999998, 0.999999964716]
for p, es, em, ep in zip(PARAMS, expected_min_s, expected_min_success, expected_safe_success):
    V = p["V"]
    s0 = minimum_s(V, 0.3)
    a0, b0 = alpha_beta(V, s0)
    a1, b1 = alpha_beta(V, s0 + 1)
    succ0 = (Decimal(1) - a0) * (Decimal(1) - b0)
    succ1 = (Decimal(1) - a1) * (Decimal(1) - b1)
    print(
        f"{p['level']:>3}: s0/s+={s0}/{s0+1}; "
        f"success={float(succ0):.12f}/{float(succ1):.12f}; "
        f"safe log2(alpha,beta)=({dec_log2(a1):.3f},{dec_log2(b1):.3f})"
    )
    if s0 != es:
        raise AssertionError(f"{p['level']} minimum field: {s0} != {es}")
    if float(succ0) < em or float(succ0) - em >= 1.1e-12:
        raise AssertionError(f"{p['level']} minimum success lower bound mismatch")
    if float(succ1) < ep or float(succ1) - ep >= 1.1e-12:
        raise AssertionError(f"{p['level']} safety success lower bound mismatch")

print("\nNORMALIZED FAST-SOLVER INDICATORS")
for eps in (0.1, 0.2, 0.3, 0.5):
    row = []
    for p in PARAMS:
        C, U = fast_indicators(p["V"], eps)
        row.append(f"{p['level']} {C:.3f}/{U:.3f} (margin {p['target']-U:.3f})")
    print(f"epsilon={eps:.1f}: " + "; ".join(row))
print("thresholds:", ", ".join(f"{epsilon_threshold(p['V'],p['target']):.6f}" for p in PARAMS))
for p, want in zip(PARAMS, [0.559284, 0.586277, 0.560058]):
    assert_close(f"{p['level']} epsilon threshold", epsilon_threshold(p["V"],p["target"]), want, 6e-7)

print("\nEPSILON-FREE SLP DIAGNOSTIC")
for p in PARAMS:
    V = p["V"]
    M2 = (V + 2) * (V + 1) // 2
    L = M2 + V * (2 * M2 - 1)
    s = ceil((log2(100) + 3 * V) / LOGQ)
    A0 = 2 * V + log2(L) + log2(s * LOGQ)
    A1 = A0 + log2(A0)
    print(f"{p['level']:>3}: s={s}, A0/A1={A0:.3f}/{A1:.3f}, margins={p['target']-A0:.3f}/{p['target']-A1:.3f}")

print("\nINDEPENDENT GIMENEZ SPECIALIZATIONS")
choices = {
    "I": ((13, 9), (8, 9)),
    "III": ((23, 13), (18, 13)),
    "V": ((7, 16), (22, 17)),
}
for p in PARAMS:
    V = p["V"]
    M = V * (V + 1) // 2
    Lcore = M + 2 * V * (M + V + 1)
    for label, (b, s) in zip(("projective", "graph-local"), choices[p["level"]]):
        eps = 2 ** (-b)
        if label == "projective":
            T = 2 * (V**4 + V * (Lcore + V**2 + V**4)) * 4**V * log2(Q**s)
        else:
            Lloc = Lcore + V**4
            T = (V**4 + V * (Lloc + V**2 + V**4)) * V * 4**V * log2(Q**s)
        direct = log2(T / (1 - 2 * eps))
        stress = direct + log2(direct)
        print(f"{p['level']:>3} {label:>11}: (2^-{b},s={s}) {direct:.3f}/{stress:.3f}")

print("\nPROJECTIVE BAD-KEY BOUNDS")
for p in PARAMS:
    V, m, O = p["V"], p["m"], p["O"]
    e = m - V
    coff = e * (O + 1) + 1
    eta_off = (
        Decimal(V) * Decimal(3) ** m * Decimal(2 * V) ** ((e + 1) * O)
        * Decimal(Q) ** (-coff)
    )
    tau = rank_failure(Q, m, V)
    eta_oil = tau ** O
    Dform = V * (V + 1) // 2 + V * O
    eta_span = gaussian_binomial(m, e, Q) * Decimal(Q) ** (-e * Dform)
    total = eta_off + eta_oil + eta_span
    print(
        f"{p['level']:>3}: coff={coff}, log2(off/oil/span/total)="
        f"{dec_log2(eta_off):.3f}/{dec_log2(eta_oil):.3f}/{dec_log2(eta_span):.3f}/{dec_log2(total):.3f}"
    )

print("\nSEEDED-MODEL TRANSFER")
for p in PARAMS:
    eps_pub = Decimal(2 * p["m"]) * Decimal(2) ** (-p["lam"])
    T = p["m"] * (ceil(p["tau1"] / 16) + ceil(p["tau2"] / 16))
    dperm = ideal_cipher_distance(T)
    dic = dperm + eps_pub
    print(
        f"{p['level']:>3}: T={T}, log2(delta_XOF/delta_IC)="
        f"{dec_log2(eps_pub):.3f}/{dec_log2(dic):.3f}"
    )

print("\nMEMORY LOWER BOUNDS (minimum epsilon=0.3 fields)")
for p in PARAMS:
    V = p["V"]
    s = minimum_s(V, 0.3)
    one = (2**V) * s * LOGQ / 8 / EIB
    arrays = V * one
    chart_log = V - s * LOGQ
    print(f"{p['level']:>3}: s={s}, one-vector={one:.12g} EiB, V arrays={arrays:.12g} EiB, log2 chart miss={chart_log:.6f}")

print("\nAll frozen numerical assertions passed.")
