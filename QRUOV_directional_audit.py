#!/usr/bin/env python3
"""Exact companion audit for the QR-UOV directional recovery paper.

The script uses only Python's standard library and performs four independent
checks:
  1. Recompute the finite cost ledger with the required 2*C0 lifting precision.
  2. Certify three sparse working-field moduli by Rabin's irreducibility test.
  3. Check the characteristic-127 start, Hensel, and Kronecker identities.
  4. Run a complete small-parameter equivalent-key recovery experiment.

The experiments validate algebraic identities; they are not full-parameter
benchmarks.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import math
import random
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

Q = 127 ** 3
LOG127 = math.log2(127)
PARAMS = [
    ("I", 52, 18, 54, 143),
    ("III", 76, 26, 78, 207),
    ("V", 102, 35, 105, 272),
]

# f(X)=X^d+sum a_i X^i over F_127.  Every nonleading coefficient is +/-1.
SPARSE_MODULI: Dict[int, Dict[int, int]] = {
    18: {4: 1, 0: 1},
    24: {5: 1, 1: -1, 0: 1},
    30: {7: 1, 0: 1},
    # These three are used only when exploring possible retuning thresholds.
    27: {12: 1, 8: 1, 0: 1},
    33: {13: -1, 0: 1},
    36: {8: 1, 1: 1, 0: 1},
}


def log2_m_from_log2_n(log2_n: float) -> float:
    """log2(n log2 n log2 log2 n), supplied log2(n)."""
    if log2_n <= 1:
        raise ValueError("stress multiplication function is used only for n>2")
    return log2_n + math.log2(log2_n) + math.log2(math.log2(log2_n))


def min_extension_degree(v: int) -> int:
    # Q^r >= 8 C^2, C=2^V.
    return math.ceil((3 + 2 * v) / math.log2(Q) - 1e-15)


def invertibility_probability(n: int, field_size: int) -> float:
    out = 1.0
    for j in range(1, n + 1):
        out *= 1.0 - field_size ** (-j)
    return out


def dense_extension_gate_envelope(d: int) -> float:
    """Old deliberately dense schoolbook multiplication/reduction envelope."""
    gmul = 2 * LOG127 * LOG127 + LOG127
    return 2 * d * d * gmul + 4 * d * d * LOG127


def schoolbook_sparse_extension_gate_envelope(d: int) -> float:
    """Schoolbook convolution plus addition-only sparse reduction."""
    lower = SPARSE_MODULI[d]
    s = len(lower)  # number of nonleading terms
    gmul = 2 * LOG127 * LOG127 + LOG127
    gadd = 4 * LOG127
    return d * d * gmul + ((d - 1) ** 2 + s * (d - 1)) * gadd


@functools.lru_cache(maxsize=None)
def karatsuba_convolution_plan(d: int) -> Tuple[float, int, int, bool]:
    """Best explicit schoolbook/Karatsuba plan for length-d convolution.

    The return value is (Boolean-gate envelope, base-field multiplications,
    base-field additions, use_Karatsuba_at_the_root).  For an uneven split
    a=ceil(d/2), b=floor(d/2), the standard three-product identity costs
    2*T(a)+T(b) and at most 4d-4 coefficient additions.  We choose the cheaper
    plan recursively in the same gate unit used by the paper.
    """
    if d <= 0:
        raise ValueError("convolution length must be positive")
    gmul = 2 * LOG127 * LOG127 + LOG127
    gadd = 4 * LOG127
    school_mul = d * d
    school_add = (d - 1) ** 2
    best = (school_mul * gmul + school_add * gadd,
            school_mul, school_add, False)
    if d >= 2:
        a = (d + 1) // 2
        b = d // 2
        pa = karatsuba_convolution_plan(a)
        pb = karatsuba_convolution_plan(b)
        kara_mul = 2 * pa[1] + pb[1]
        kara_add = 2 * pa[2] + pb[2] + 4 * d - 4
        candidate = (kara_mul * gmul + kara_add * gadd,
                     kara_mul, kara_add, True)
        if candidate[0] < best[0]:
            best = candidate
    return best


def sparse_extension_gate_envelope(d: int) -> float:
    """Karatsuba convolution plus addition-only reduction by the modulus."""
    lower = SPARSE_MODULI[d]
    s = len(lower)
    gadd = 4 * LOG127
    convolution, _, _, _ = karatsuba_convolution_plan(d)
    return convolution + s * (d - 1) * gadd


@dataclass
class CostRow:
    level: str
    V: int
    O: int
    m: int
    target: int
    r: int
    d: int
    sparse_modulus: str
    restart_factor: float
    separated_factor_stress: float
    corrected_dense_field_stress: float
    corrected_schoolbook_sparse_field_stress: float
    corrected_sparse_field_stress: float
    sparse_margin: float
    extension_saving_bits: float
    karatsuba_saving_bits: float
    one_flat_array_memory_log2_bytes: float


def modulus_string(d: int) -> str:
    terms = [f"X^{d}"]
    for e in sorted(SPARSE_MODULI[d], reverse=True):
        a = SPARSE_MODULI[d][e]
        mon = "1" if e == 0 else ("X" if e == 1 else f"X^{e}")
        terms.append(("+" if a == 1 else "-") + mon)
    return "".join(terms)


def cost_row(level: str, V: int, O: int, m: int, target: int) -> CostRow:
    log_c = float(V)
    log_c0 = V - 1 + math.log2(V + 2)
    r = min_extension_degree(V)
    d = 3 * r
    if d not in SPARSE_MODULI:
        raise KeyError(f"no audited sparse modulus for degree {d}")

    p_dir = 1.0 - V / Q
    p_mix = invertibility_probability(V, Q)
    p_sep = 1.0 - 2.0 ** (2 * V - r * math.log2(Q))
    if p_sep <= 0:
        raise AssertionError("working extension does not meet separator bound")
    restart = 1.0 / (p_dir * p_mix * p_sep)

    log_outer = math.log2(V ** 3 + V ** 2)
    dense_ext = dense_extension_gate_envelope(d)
    schoolbook_sparse_ext = schoolbook_sparse_extension_gate_envelope(d)
    sparse_ext = sparse_extension_gate_envelope(d)
    log_restart = math.log2(restart)

    # Legacy separated-axis stress, retained only as a comparison.
    separated = (log_outer + math.log2(dense_ext) + log_restart
                 + log2_m_from_log2_n(log_c)
                 + log2_m_from_log2_n(log_c0))

    # The global curve has degree C0 and rational reconstruction requires
    # precision 2*C0.  At k=2*C0, a flattened input has length <4*C*C0;
    # a full multiply-and-reduce schedule is conservatively charged by
    # 3*M(8*C*C0).  This includes the correction omitted in the previous draft.
    log_schedule_n = 3 + log_c + log_c0
    dense = (log_outer + math.log2(dense_ext) + log_restart + math.log2(3)
             + log2_m_from_log2_n(log_schedule_n))
    schoolbook_sparse = (log_outer + math.log2(schoolbook_sparse_ext)
                         + log_restart + math.log2(3)
                         + log2_m_from_log2_n(log_schedule_n))
    sparse = (log_outer + math.log2(sparse_ext) + log_restart + math.log2(3)
              + log2_m_from_log2_n(log_schedule_n))

    # One coefficient array has C*(2*C0) elements of F_{127^d}.
    log_mem_bytes = log_c + 1 + log_c0 + math.log2(d * LOG127) - 3

    return CostRow(
        level=level, V=V, O=O, m=m, target=target, r=r, d=d,
        sparse_modulus=modulus_string(d), restart_factor=restart,
        separated_factor_stress=separated,
        corrected_dense_field_stress=dense,
        corrected_schoolbook_sparse_field_stress=schoolbook_sparse,
        corrected_sparse_field_stress=sparse,
        sparse_margin=target - sparse,
        extension_saving_bits=math.log2(dense_ext / sparse_ext),
        karatsuba_saving_bits=math.log2(schoolbook_sparse_ext / sparse_ext),
        one_flat_array_memory_log2_bytes=log_mem_bytes,
    )

# ---------------------------------------------------------------------------
# Elementary finite-field and power-series routines
# ---------------------------------------------------------------------------


def inv_mod(a: int, p: int) -> int:
    a %= p
    if a == 0:
        raise ZeroDivisionError
    return pow(a, p - 2, p)


def mat_inv_mod(a: Sequence[Sequence[int]], p: int) -> List[List[int]]:
    n = len(a)
    aug = [list(row[i] % p for i in range(n)) + [1 if i == j else 0 for j in range(n)]
           for i, row in enumerate(a)]
    for c in range(n):
        pivot = next((r for r in range(c, n) if aug[r][c] % p), None)
        if pivot is None:
            raise ZeroDivisionError("singular matrix")
        aug[c], aug[pivot] = aug[pivot], aug[c]
        s = inv_mod(aug[c][c], p)
        aug[c] = [(s * x) % p for x in aug[c]]
        for r in range(n):
            if r == c:
                continue
            f = aug[r][c] % p
            if f:
                aug[r] = [(x - f * y) % p for x, y in zip(aug[r], aug[c])]
    return [row[n:] for row in aug]


def det_mod(a: Sequence[Sequence[int]], p: int) -> int:
    m = [list(x % p for x in row) for row in a]
    n = len(m)
    det = 1
    for c in range(n):
        pivot = next((r for r in range(c, n) if m[r][c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            m[c], m[pivot] = m[pivot], m[c]
            det = -det
        piv = m[c][c] % p
        det = det * piv % p
        ip = inv_mod(piv, p)
        for r in range(c + 1, n):
            f = m[r][c] * ip % p
            for j in range(c, n):
                m[r][j] = (m[r][j] - f * m[c][j]) % p
    return det % p


def poly_from_roots(roots: Sequence[int], p: int) -> List[int]:
    # Increasing coefficient order.
    out = [1]
    for a in roots:
        nxt = [0] * (len(out) + 1)
        for i, c in enumerate(out):
            nxt[i] = (nxt[i] - a * c) % p
            nxt[i + 1] = (nxt[i + 1] + c) % p
        out = nxt
    return out


def series_add(a: Sequence[int], b: Sequence[int], p: int, k: int) -> List[int]:
    return [((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(k)]


def series_sub(a: Sequence[int], b: Sequence[int], p: int, k: int) -> List[int]:
    return [((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(k)]


def series_mul(a: Sequence[int], b: Sequence[int], p: int, k: int) -> List[int]:
    out = [0] * k
    for i, x in enumerate(a[:k]):
        if x:
            for j, y in enumerate(b[:k - i]):
                out[i + j] = (out[i + j] + x * y) % p
    return out


def series_scale(a: Sequence[int], c: int, p: int, k: int) -> List[int]:
    return [(c * (a[i] if i < len(a) else 0)) % p for i in range(k)]


def series_matrix_mul(A, B, p: int, k: int):
    n, mid, m = len(A), len(B), len(B[0])
    out = [[[0] * k for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for r in range(mid):
            for j in range(m):
                out[i][j] = series_add(out[i][j], series_mul(A[i][r], B[r][j], p, k), p, k)
    return out


def series_matrix_inverse(A, p: int, k: int):
    n = len(A)
    A0 = [[A[i][j][0] % p for j in range(n)] for i in range(n)]
    B0 = mat_inv_mod(A0, p)
    B = [[[B0[i][j]] for j in range(n)] for i in range(n)]
    precision = 1
    while precision < k:
        newk = min(2 * precision, k)
        Bt = [[[x[t] if t < len(x) else 0 for t in range(newk)] for x in row] for row in B]
        At = [[[x[t] if t < len(x) else 0 for t in range(newk)] for x in row] for row in A]
        AB = series_matrix_mul(At, Bt, p, newk)
        twoI_minus = [[[0] * newk for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                twoI_minus[i][j][0] = ((2 if i == j else 0) - AB[i][j][0]) % p
                for t in range(1, newk):
                    twoI_minus[i][j][t] = (-AB[i][j][t]) % p
        B = series_matrix_mul(Bt, twoI_minus, p, newk)
        precision = newk
    return B


@dataclass
class Quad:
    # f(x)=c+b^T x+sum_{a,b} A[a][b] x_a x_b, no symmetry required.
    A: List[List[int]]
    b: List[int]
    c: int


def quad_eval_series(q: Quad, X: Sequence[Sequence[int]], p: int, k: int) -> List[int]:
    v = len(X)
    out = [q.c % p] + [0] * (k - 1)
    for i in range(v):
        out = series_add(out, series_scale(X[i], q.b[i], p, k), p, k)
        for j in range(v):
            if q.A[i][j] % p:
                out = series_add(out, series_scale(series_mul(X[i], X[j], p, k), q.A[i][j], p, k), p, k)
    return out


def quad_jac_series(q: Quad, X: Sequence[Sequence[int]], p: int, k: int) -> List[List[int]]:
    v = len(X)
    row = []
    for j in range(v):
        s = [q.b[j] % p] + [0] * (k - 1)
        for r in range(v):
            coeff = (q.A[j][r] + q.A[r][j]) % p
            if coeff:
                s = series_add(s, series_scale(X[r], coeff, p, k), p, k)
        row.append(s)
    return row


def affine_linear_product_quad(a: int, b: int, V: int, p: int) -> Quad:
    # L_a=X1+aX2+...+a^(V-1)XV+a^V.
    la = [pow(a, j, p) for j in range(V)]
    lb = [pow(b, j, p) for j in range(V)]
    ca, cb = pow(a, V, p), pow(b, V, p)
    A = [[la[i] * lb[j] % p for j in range(V)] for i in range(V)]
    lin = [(ca * lb[i] + cb * la[i]) % p for i in range(V)]
    return Quad(A=A, b=lin, c=ca * cb % p)


def random_quad(V: int, p: int, rng: random.Random) -> Quad:
    A = [[rng.randrange(p) for _ in range(V)] for _ in range(V)]
    b = [rng.randrange(p) for _ in range(V)]
    return Quad(A=A, b=b, c=rng.randrange(p))


def start_root(selected_nodes: Sequence[int], p: int) -> List[int]:
    coeff = poly_from_roots(selected_nodes, p)  # c0,...,cV with cV=1
    V = len(selected_nodes)
    # Polynomial is T^V + X_V T^(V-1)+...+X_1.
    return [coeff[j] % p for j in range(V)]


def homotopy_eval(g: Sequence[Quad], f: Sequence[Quad], X, p: int, k: int):
    t = [0, 1] + [0] * max(0, k - 2)
    one_minus_t = [1, -1 % p] + [0] * max(0, k - 2)
    out = []
    for gi, fi in zip(g, f):
        gv = quad_eval_series(gi, X, p, k)
        fv = quad_eval_series(fi, X, p, k)
        out.append(series_add(series_mul(one_minus_t, gv, p, k), series_mul(t, fv, p, k), p, k))
    return out


def homotopy_jac(g: Sequence[Quad], f: Sequence[Quad], X, p: int, k: int):
    t = [0, 1] + [0] * max(0, k - 2)
    one_minus_t = [1, -1 % p] + [0] * max(0, k - 2)
    out = []
    for gi, fi in zip(g, f):
        gj = quad_jac_series(gi, X, p, k)
        fj = quad_jac_series(fi, X, p, k)
        out.append([series_add(series_mul(one_minus_t, a, p, k), series_mul(t, b, p, k), p, k)
                    for a, b in zip(gj, fj)])
    return out


def lift_coefficientwise(g, f, x0, p: int, precision: int):
    V = len(x0)
    X = [[x0[i]] + [0] * (precision - 1) for i in range(V)]
    J0 = [[homotopy_jac(g, f, [[x0[j]] for j in range(V)], p, 1)[i][j][0]
           for j in range(V)] for i in range(V)]
    J0inv = mat_inv_mod(J0, p)
    for r in range(1, precision):
        H = homotopy_eval(g, f, X, p, r + 1)
        residual = [H[i][r] % p for i in range(V)]
        delta = [(-sum(J0inv[i][j] * residual[j] for j in range(V))) % p for i in range(V)]
        for i in range(V):
            X[i][r] = delta[i]
    return X


def lift_newton(g, f, x0, p: int, precision: int):
    V = len(x0)
    X = [[x0[i]] for i in range(V)]
    cur = 1
    while cur < precision:
        newk = min(2 * cur, precision)
        X = [xi + [0] * (newk - len(xi)) for xi in X]
        H = homotopy_eval(g, f, X, p, newk)
        J = homotopy_jac(g, f, X, p, newk)
        Jinv = series_matrix_inverse(J, p, newk)
        Hcol = [[h] for h in H]
        correction = series_matrix_mul(Jinv, Hcol, p, newk)
        X = [series_sub(X[i], correction[i][0], p, newk) for i in range(V)]
        cur = newk
    return X


# ---------------------------------------------------------------------------
# Bivariate multiplication / flattening audit
# ---------------------------------------------------------------------------


def tconv(a: Sequence[int], b: Sequence[int], p: int, k: int) -> List[int]:
    return series_mul(a, b, p, k)


def direct_u_product(A, B, p: int, k: int):
    out = [[0] * k for _ in range(len(A) + len(B) - 1)]
    for i, ai in enumerate(A):
        for j, bj in enumerate(B):
            out[i + j] = series_add(out[i + j], tconv(ai, bj, p, k), p, k)
    return out


def flatten_u_product(A, B, p: int, k: int):
    C = max(len(A), len(B))
    radix = 2 * C - 1
    def encode(P):
        n = radix * k
        out = [0] * n
        for u, coeff in enumerate(P):
            for t, x in enumerate(coeff[:k]):
                out[u + radix * t] = x % p
        return out
    ea, eb = encode(A), encode(B)
    raw = [0] * (len(ea) + len(eb) - 1)
    for i, x in enumerate(ea):
        if x:
            for j, y in enumerate(eb):
                if y:
                    raw[i + j] = (raw[i + j] + x * y) % p
    out = [[0] * k for _ in range(len(A) + len(B) - 1)]
    for e, x in enumerate(raw):
        if not x:
            continue
        t, u = divmod(e, radix)
        if t < k and u < len(out):
            out[u][t] = (out[u][t] + x) % p
    return out


def reduce_monic_u(P, Qpoly, p: int, k: int):
    # Qpoly has U degree C and leading coefficient 1.
    C = len(Qpoly) - 1
    R = [row[:] for row in P] + [[0] * k for _ in range(max(0, 2 * C - 1 - len(P)))]
    for u in range(len(R) - 1, C - 1, -1):
        high = R[u][:]
        if any(high):
            for j in range(C):
                prod = tconv(high, Qpoly[j], p, k)
                R[u - C + j] = series_sub(R[u - C + j], prod, p, k)
            R[u] = [0] * k
    return R[:C]


def random_bivariate_poly(C: int, k: int, p: int, rng: random.Random):
    return [[rng.randrange(p) for _ in range(k)] for _ in range(C)]



# ---------------------------------------------------------------------------
# Polynomial arithmetic and Rabin irreducibility certificates
# ---------------------------------------------------------------------------


def poly_trim(a: Sequence[int], p: int) -> List[int]:
    out = [x % p for x in a]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_add_mod(a: Sequence[int], b: Sequence[int], p: int) -> List[int]:
    n = max(len(a), len(b))
    return poly_trim([((a[i] if i < len(a) else 0) +
                       (b[i] if i < len(b) else 0)) % p for i in range(n)], p)


def poly_sub_mod(a: Sequence[int], b: Sequence[int], p: int) -> List[int]:
    n = max(len(a), len(b))
    return poly_trim([((a[i] if i < len(a) else 0) -
                       (b[i] if i < len(b) else 0)) % p for i in range(n)], p)


def poly_mul_plain(a: Sequence[int], b: Sequence[int], p: int) -> List[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return poly_trim(out, p)


def _poly_mul_planned_equal(a: Sequence[int], b: Sequence[int], p: int) -> List[int]:
    """Execute the recursive plan certified by karatsuba_convolution_plan.

    Inputs have equal positive length.  Zero padding is kept internally so the
    operation follows the same uneven-split recurrence used in the gate count.
    """
    if len(a) != len(b) or not a:
        raise ValueError("planned multiplication expects equal positive lengths")
    n = len(a)
    if not karatsuba_convolution_plan(n)[3]:
        out = [0] * (2 * n - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
        return out

    h = (n + 1) // 2
    t = n // 2
    a0 = list(a[:h])
    b0 = list(b[:h])
    a1 = list(a[h:])
    b1 = list(b[h:])
    z0 = _poly_mul_planned_equal(a0, b0, p)
    z2 = _poly_mul_planned_equal(a1, b1, p)
    sa = [(a0[i] + (a1[i] if i < t else 0)) % p for i in range(h)]
    sb = [(b0[i] + (b1[i] if i < t else 0)) % p for i in range(h)]
    z1 = _poly_mul_planned_equal(sa, sb, p)
    for i in range(len(z0)):
        z1[i] = (z1[i] - z0[i]) % p
    for i in range(len(z2)):
        z1[i] = (z1[i] - z2[i]) % p

    out = [0] * (2 * n - 1)
    for i, x in enumerate(z0):
        out[i] = (out[i] + x) % p
    for i, x in enumerate(z1):
        if i + h < len(out):
            out[i + h] = (out[i + h] + x) % p
    for i, x in enumerate(z2):
        if i + 2 * h < len(out):
            out[i + 2 * h] = (out[i + 2 * h] + x) % p
    return out


def poly_mul_planned(a: Sequence[int], b: Sequence[int], p: int) -> List[int]:
    """Multiply arbitrary nonempty coefficient vectors with the audited plan."""
    if not a or not b:
        return [0]
    n = max(len(a), len(b))
    aa = list(a) + [0] * (n - len(a))
    bb = list(b) + [0] * (n - len(b))
    out = _poly_mul_planned_equal(aa, bb, p)[:len(a) + len(b) - 1]
    return poly_trim(out, p)


def poly_divmod_mod(a: Sequence[int], b: Sequence[int], p: int):
    r = poly_trim(a, p)
    b = poly_trim(b, p)
    if b == [0]:
        raise ZeroDivisionError
    q = [0] * max(1, len(r) - len(b) + 1)
    ib = inv_mod(b[-1], p)
    while len(r) >= len(b) and r != [0]:
        k = len(r) - len(b)
        c = r[-1] * ib % p
        q[k] = c
        for j in range(len(b)):
            r[k + j] = (r[k + j] - c * b[j]) % p
        r = poly_trim(r, p)
    return poly_trim(q, p), r


def poly_mod(a: Sequence[int], f: Sequence[int], p: int) -> List[int]:
    return poly_divmod_mod(a, f, p)[1]


def poly_mul_mod(a: Sequence[int], b: Sequence[int], f: Sequence[int], p: int) -> List[int]:
    return poly_mod(poly_mul_plain(a, b, p), f, p)


def poly_pow_mod(a: Sequence[int], e: int, f: Sequence[int], p: int) -> List[int]:
    out = [1]
    base = poly_mod(a, f, p)
    while e:
        if e & 1:
            out = poly_mul_mod(out, base, f, p)
        base = poly_mul_mod(base, base, f, p)
        e >>= 1
    return out


def poly_gcd_mod(a: Sequence[int], b: Sequence[int], p: int) -> List[int]:
    a, b = poly_trim(a, p), poly_trim(b, p)
    while b != [0]:
        _, r = poly_divmod_mod(a, b, p)
        a, b = b, r
    inv = inv_mod(a[-1], p)
    return [(inv * x) % p for x in a]


def prime_divisors(n: int) -> List[int]:
    out = []
    q = 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out


def sparse_modulus_coeffs(d: int, p: int = 127) -> List[int]:
    f = [0] * (d + 1)
    f[d] = 1
    for e, a in SPARSE_MODULI[d].items():
        f[e] = a % p
    return f


def rabin_irreducibility_certificate(d: int, p: int = 127) -> dict:
    f = sparse_modulus_coeffs(d, p)
    x = [0, 1]
    frob = {0: x}
    h = x
    for k in range(1, d + 1):
        h = poly_pow_mod(h, p, f, p)
        frob[k] = h
    final_ok = poly_sub_mod(frob[d], x, p) == [0]
    gcd_degrees = {}
    for ell in prime_divisors(d):
        g = poly_gcd_mod(f, poly_sub_mod(frob[d // ell], x, p), p)
        gcd_degrees[str(ell)] = len(g) - 1
    irreducible = final_ok and all(v == 0 for v in gcd_degrees.values())
    if not irreducible:
        raise AssertionError(f"sparse modulus degree {d} failed Rabin test")
    return {
        "degree": d,
        "polynomial": modulus_string(d),
        "x_qd_equals_x": final_ok,
        "gcd_degrees_for_prime_divisors": gcd_degrees,
        "irreducible": irreducible,
    }


# ---------------------------------------------------------------------------
# Complete toy equivalent-key recovery
# ---------------------------------------------------------------------------


def mat_vec(A, x, p):
    return [sum(a * b for a, b in zip(row, x)) % p for row in A]


def vec_dot(a, b, p):
    return sum(x * y for x, y in zip(a, b)) % p


def mat_mul(A, B, p):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % p
             for j in range(len(B[0]))] for i in range(len(A))]


def mat_transpose(A):
    return [list(x) for x in zip(*A)]


def mat_add(A, B, p):
    return [[(x + y) % p for x, y in zip(r, s)] for r, s in zip(A, B)]


def mat_sub(A, B, p):
    return [[(x - y) % p for x, y in zip(r, s)] for r, s in zip(A, B)]


def solve_full_column_rank(M, b, p):
    rows, cols = len(M), len(M[0])
    aug = [[M[i][j] % p for j in range(cols)] + [b[i] % p] for i in range(rows)]
    pivots = []
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if aug[i][c]), None)
        if piv is None:
            continue
        aug[r], aug[piv] = aug[piv], aug[r]
        s = inv_mod(aug[r][c], p)
        aug[r] = [(s * z) % p for z in aug[r]]
        for i in range(rows):
            if i != r and aug[i][c]:
                t = aug[i][c]
                aug[i] = [(u - t * v) % p for u, v in zip(aug[i], aug[r])]
        pivots.append(c)
        r += 1
    for i in range(r, rows):
        if all(aug[i][j] == 0 for j in range(cols)) and aug[i][-1] != 0:
            return None
    if len(pivots) != cols:
        return None
    x = [0] * cols
    for i, c in enumerate(pivots):
        x[c] = aug[i][-1]
    return x


def random_symmetric(n, p, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            A[i][j] = A[j][i] = rng.randrange(p)
    return A


def toy_public_key(V, O, m, p, rng):
    G = [[rng.randrange(p) for _ in range(O)] for _ in range(V)]
    As, Es, Cs = [], [], []
    for _ in range(m):
        A = random_symmetric(V, p, rng)
        E = [[rng.randrange(p) for _ in range(O)] for _ in range(V)]
        GT, ET = mat_transpose(G), mat_transpose(E)
        C = mat_sub(mat_add(mat_mul(GT, E, p), mat_mul(ET, G, p), p),
                    mat_mul(mat_mul(GT, A, p), G, p), p)
        As.append(A); Es.append(E); Cs.append(C)
    return G, As, Es, Cs


def directional_value(A, E, C, x, c, p):
    return (vec_dot(x, mat_vec(A, x, p), p)
            - 2 * vec_dot(x, mat_vec(E, c, p), p)
            + vec_dot(c, mat_vec(C, c, p), p)) % p


def directional_jacobian(As, Es, x, c, p):
    return [[2 * (u - v) % p for u, v in zip(mat_vec(A, x, p), mat_vec(E, c, p))]
            for A, E in zip(As, Es)]


def recover_graph_from_root(As, Es, Cs, g, c, O, p):
    V = len(g)
    cols = []
    for j in range(O):
        d = [1 if t == j else 0 for t in range(O)]
        M, rhs = [], []
        for A, E, C in zip(As, Es, Cs):
            M.append([(u - v) % p for u, v in zip(mat_vec(A, g, p), mat_vec(E, c, p))])
            rhs.append((vec_dot(g, mat_vec(E, d, p), p)
                        - vec_dot(c, mat_vec(C, d, p), p)) % p)
        y = solve_full_column_rank(M, rhs, p)
        if y is None:
            return None
        cols.append(y)
    return [[cols[j][i] for j in range(O)] for i in range(V)]


def verify_graph(G, As, Es, Cs, p):
    GT = mat_transpose(G)
    for A, E, C in zip(As, Es, Cs):
        expected = mat_sub(mat_add(mat_mul(GT, E, p), mat_mul(mat_transpose(E), G, p), p),
                           mat_mul(mat_mul(GT, A, p), G, p), p)
        if expected != C:
            return False
    return True


def run_toy_equivalent_key_experiment(seed: int, trials: int):
    rng = random.Random(seed ^ 0x5152554F56)
    p, V, O, m = 7, 3, 2, 5
    successful_directions = 0
    recovered = 0
    false_candidates = 0
    roots_examined = 0
    trials_without_verified_key = 0
    verified_keys_from_nonplanted_roots = 0
    for _ in range(trials):
        G, As, Es, Cs = toy_public_key(V, O, m, p, rng)
        c = [rng.randrange(p) for _ in range(O)]
        if all(z == 0 for z in c):
            c[0] = 1
        planted = mat_vec(G, c, p)
        roots = []
        for x in itertools.product(range(p), repeat=V):
            if all(directional_value(A, E, C, x, c, p) == 0
                   for A, E, C in zip(As, Es, Cs)):
                roots.append(list(x))
        if planted not in roots:
            raise AssertionError("planted directional root disappeared")
        roots_examined += len(roots)
        accepted = []
        for root in roots:
            J = directional_jacobian(As, Es, root, c, p)
            candidate = recover_graph_from_root(As, Es, Cs, root, c, O, p)
            if candidate is not None and verify_graph(candidate, As, Es, Cs, p):
                accepted.append(candidate)
            elif candidate is not None:
                false_candidates += 1
        planted_rank_candidate = recover_graph_from_root(As, Es, Cs, planted, c, O, p)
        if planted_rank_candidate is not None:
            successful_directions += 1
            if planted_rank_candidate != G:
                raise AssertionError("completion at planted root returned wrong graph")
            if not any(A == G for A in accepted):
                raise AssertionError("exact verification rejected the true graph")
            recovered += 1
            verified_keys_from_nonplanted_roots += sum(A != G for A in accepted)
        else:
            if accepted:
                verified_keys_from_nonplanted_roots += len(accepted)
            else:
                trials_without_verified_key += 1
    if recovered == 0:
        raise AssertionError("toy experiment recovered no equivalent key")
    return {
        "field": p, "V": V, "O": O, "m": m, "trials": trials,
        "full_rank_planted_directions": successful_directions,
        "verified_equivalent_keys_recovered": recovered,
        "all_equation_roots_examined": roots_examined,
        "nonverifying_full_rank_candidates": false_candidates,
        "trials_without_verified_key": trials_without_verified_key,
        "verified_keys_from_nonplanted_roots": verified_keys_from_nonplanted_roots,
    }


# ---------------------------------------------------------------------------
# Test driver
# ---------------------------------------------------------------------------


def run_exact_tests(seed: int, trials: int):
    rng = random.Random(seed)
    start_checks = 0
    lifted_paths = 0
    flatten_checks = 0
    karatsuba_checks = 0

    for V in (2, 3, 4):
        p = 127
        nodes = list(range(1, 2 * V + 1))
        g = [affine_linear_product_quad(nodes[2 * i], nodes[2 * i + 1], V, p) for i in range(V)]
        for choice in itertools.product((0, 1), repeat=V):
            selected = [nodes[2 * i + choice[i]] for i in range(V)]
            x0 = start_root(selected, p)
            vals = [quad_eval_series(q, [[x] for x in x0], p, 1)[0] for q in g]
            if any(vals):
                raise AssertionError("start root does not satisfy product system")
            J = [[quad_jac_series(q, [[x] for x in x0], p, 1)[j][0]
                  for j in range(V)] for q in g]
            if det_mod(J, p) == 0:
                raise AssertionError("start root is singular")
            start_checks += 1

        for _ in range(trials):
            f = [random_quad(V, p, rng) for _ in range(V)]
            choice = [rng.randrange(2) for _ in range(V)]
            selected = [nodes[2 * i + choice[i]] for i in range(V)]
            x0 = start_root(selected, p)
            precision = 10
            a = lift_coefficientwise(g, f, x0, p, precision)
            b = lift_newton(g, f, x0, p, precision)
            if a != b:
                raise AssertionError("coefficientwise and Newton lifts disagree")
            if any(any(c % p for c in h) for h in homotopy_eval(g, f, b, p, precision)):
                raise AssertionError("lifted path has nonzero homotopy residual")
            lifted_paths += 1

    for C in (3, 4, 5):
        for k in (3, 5, 7):
            p = 127
            for _ in range(trials):
                A = random_bivariate_poly(C, k, p, rng)
                B = random_bivariate_poly(C, k, p, rng)
                Qpoly = random_bivariate_poly(C, k, p, rng) + [[1] + [0] * (k - 1)]
                direct = direct_u_product(A, B, p, k)
                flat = flatten_u_product(A, B, p, k)
                if direct != flat:
                    raise AssertionError("Kronecker flattening disagrees with direct product")
                if reduce_monic_u(direct, Qpoly, p, k) != reduce_monic_u(flat, Qpoly, p, k):
                    raise AssertionError("reduced flattened product disagrees")
                flatten_checks += 1

    # Check the exact recursive convolution schedule used by every listed
    # auxiliary field against elementary schoolbook multiplication.
    p = 127
    for d in sorted(SPARSE_MODULI):
        for _ in range(max(4, trials // 5)):
            A = [rng.randrange(p) for _ in range(d)]
            B = [rng.randrange(p) for _ in range(d)]
            if poly_mul_planned(A, B, p) != poly_mul_plain(A, B, p):
                raise AssertionError(f"Karatsuba plan disagrees at degree {d}")
            karatsuba_checks += 1

    return {
        "vandermonde_start_roots_checked": start_checks,
        "characteristic_127_lifted_paths_checked": lifted_paths,
        "flattened_bivariate_products_checked": flatten_checks,
        "karatsuba_convolutions_checked": karatsuba_checks,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--trials", type=int, default=40)
    ap.add_argument("--toy-trials", type=int, default=200)
    ap.add_argument("--json", type=Path, default=Path("QRUOV_verdict_change_audit.json"))
    args = ap.parse_args()

    rows = [cost_row(*p) for p in PARAMS]
    moduli = [rabin_irreducibility_certificate(d) for d in sorted(SPARSE_MODULI)]
    identities = run_exact_tests(args.seed, args.trials)
    toy = run_toy_equivalent_key_experiment(args.seed, args.toy_trials)

    expected = {
        "I": (156.270, 154.826),
        "III": (207.891, 206.177),
        "V": (262.734, 260.857),
    }
    for row in rows:
        dense, sparse = expected[row.level]
        if abs(row.corrected_dense_field_stress - dense) > 5e-4:
            raise AssertionError(f"dense headline drift at Level {row.level}")
        if abs(row.corrected_sparse_field_stress - sparse) > 5e-4:
            raise AssertionError(f"sparse headline drift at Level {row.level}")
    if args.seed == 20260728 and args.trials == 40 and args.toy_trials == 200:
        expected_ids = {
            "vandermonde_start_roots_checked": 28,
            "characteristic_127_lifted_paths_checked": 120,
            "flattened_bivariate_products_checked": 360,
            "karatsuba_convolutions_checked": 48,
        }
        if identities != expected_ids:
            raise AssertionError("exact-identity regression count drift")
        if (toy["verified_equivalent_keys_recovered"],
            toy["trials_without_verified_key"],
            toy["verified_keys_from_nonplanted_roots"]) != (199, 1, 0):
            raise AssertionError("end-to-end toy regression drift")

    payload = {
        "description": "Corrected finite-precision, sparse-field, and explicit Karatsuba ledger for the QR-UOV directional attack.",
        "cost_rows": [asdict(r) for r in rows],
        "sparse_modulus_certificates": moduli,
        "karatsuba_convolution_plans": [
            {
                "degree": d,
                "base_field_multiplications": karatsuba_convolution_plan(d)[1],
                "base_field_additions": karatsuba_convolution_plan(d)[2],
                "schoolbook_gate_envelope": schoolbook_sparse_extension_gate_envelope(d),
                "karatsuba_gate_envelope": sparse_extension_gate_envelope(d),
                "saving_bits": math.log2(
                    schoolbook_sparse_extension_gate_envelope(d)
                    / sparse_extension_gate_envelope(d)
                ),
            }
            for d in sorted(SPARSE_MODULI)
        ],
        "exact_identity_checks": identities,
        "end_to_end_toy_equivalent_key_recovery": toy,
    }
    args.json.write_text(json.dumps(payload, indent=2) + "\n")

    print("QR-UOV directional equivalent-key audit")
    print("=" * 62)
    for r in rows:
        print(
            f"Level {r.level}: r={r.r}, d={r.d}, restart={r.restart_factor:.9f}, "
            f"separated={r.separated_factor_stress:.3f}, "
            f"corrected-dense={r.corrected_dense_field_stress:.3f}, "
            f"schoolbook-sparse={r.corrected_schoolbook_sparse_field_stress:.3f}, "
            f"Karatsuba-sparse={r.corrected_sparse_field_stress:.3f} "
            f"(target margin {r.sparse_margin:+.3f}), "
            f"total field saving={r.extension_saving_bits:.3f} bits, "
            f"Karatsuba saving={r.karatsuba_saving_bits:.3f} bits"
        )
    for c in moduli:
        print("Rabin certificate:", c["polynomial"], c["gcd_degrees_for_prime_divisors"])
    print("Exact identities:", identities)
    print("End-to-end toy recovery:", toy)
    print("Wrote", args.json)


if __name__ == "__main__":
    main()
