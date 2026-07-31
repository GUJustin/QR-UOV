#!/usr/bin/env python3
"""Reduced-parameter end-to-end QR-UOV local-separator forgery.

This program instantiates the candidate attack pipeline rather than merely
checking its component lemmas.  It implements:

* compact degree-3 QR-UOV key generation over K = F_{q^3};
* a SHAKE256-based reduced signing/verification interface;
* target normalization and the anchored Chevalley--Warning construction;
* the a = m-4 affine projective core;
* Vandermonde product-start homotopy with coefficientwise branch lifting;
* terminal local-separator characteristic-polynomial recombination;
* rational reconstruction and valuation-safe target specialization;
* exact one-sketch Frobenius descent;
* one on-demand coordinate rerun; and
* affine rescaling followed by unchanged public verification.

The default profile q=13, ell=3, V=3, O=2 has n=15, m=6 and a=2.
It is deliberately tiny, but it follows the same algorithmic path claimed in
the paper.  It is not a performance extrapolation to the NIST parameters.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import random
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

_HELPER = Path(__file__).resolve().with_name("QRUOV_cw_helper.py")
spec = importlib.util.spec_from_file_location("qruov_cw_helper", _HELPER)
h = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = h
assert spec.loader is not None
spec.loader.exec_module(h)

Vector = List[int]
Matrix = List[List[int]]
Elt3 = Tuple[int, int, int]
KMatrix = List[List[Elt3]]
KVector = List[Elt3]


# ---------------------------------------------------------------------------
# K = F_{q^3} and reduced compact QR-UOV
# ---------------------------------------------------------------------------


class GFq3:
    """F_q[z]/(z^3-z-1), for primes q where the cubic is irreducible."""

    def __init__(self, q: int):
        self.q = q
        if any((a**3 - a - 1) % q == 0 for a in range(q)):
            raise ValueError(f"z^3-z-1 is reducible over F_{q}")

    def z(self) -> Elt3:
        return (0, 0, 0)

    def o(self) -> Elt3:
        return (1, 0, 0)

    def add(self, a: Elt3, b: Elt3) -> Elt3:
        q = self.q
        return ((a[0] + b[0]) % q, (a[1] + b[1]) % q, (a[2] + b[2]) % q)

    def neg(self, a: Elt3) -> Elt3:
        q = self.q
        return ((-a[0]) % q, (-a[1]) % q, (-a[2]) % q)

    def sub(self, a: Elt3, b: Elt3) -> Elt3:
        return self.add(a, self.neg(b))

    def mul(self, a: Elt3, b: Elt3) -> Elt3:
        q = self.q
        c = [0] * 5
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i + j] = (c[i + j] + x * y) % q
        # z^3=z+1 and z^4=z^2+z.
        for d in range(4, 2, -1):
            x = c[d] % q
            if x:
                c[d] = 0
                c[d - 2] = (c[d - 2] + x) % q
                c[d - 3] = (c[d - 3] + x) % q
        return (c[0] % q, c[1] % q, c[2] % q)

    def pow(self, a: Elt3, e: int) -> Elt3:
        out = self.o()
        base = a
        while e:
            if e & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            e >>= 1
        return out

    def inv(self, a: Elt3) -> Elt3:
        if a == self.z():
            raise ZeroDivisionError
        return self.pow(a, self.q**3 - 2)

    def rand(self, rng: random.Random) -> Elt3:
        return (rng.randrange(self.q), rng.randrange(self.q), rng.randrange(self.q))

    def trace(self, a: Elt3) -> int:
        aq = self.pow(a, self.q)
        aq2 = self.pow(aq, self.q)
        return self.add(self.add(a, aq), aq2)[0] % self.q


def kdot(F: GFq3, x: KVector, y: KVector) -> Elt3:
    s = F.z()
    for a, b in zip(x, y):
        s = F.add(s, F.mul(a, b))
    return s


def kmat_vec(F: GFq3, A: KMatrix, x: KVector) -> KVector:
    return [kdot(F, row, x) for row in A]


def kmat_mul(F: GFq3, A: KMatrix, B: KMatrix) -> KMatrix:
    Bt = [list(col) for col in zip(*B)]
    return [[kdot(F, row, col) for col in Bt] for row in A]


def ktranspose(A: KMatrix) -> KMatrix:
    return [list(col) for col in zip(*A)]


def kquad(F: GFq3, A: KMatrix, x: KVector) -> Elt3:
    return kdot(F, x, kmat_vec(F, A, x))


def base_to_k(x: Vector) -> KVector:
    assert len(x) % 3 == 0
    return [tuple(x[3 * i : 3 * i + 3]) for i in range(len(x) // 3)]  # type: ignore[list-item]


def k_to_base(x: KVector) -> Vector:
    return [c for a in x for c in a]


def random_ksymmetric(F: GFq3, n: int, rng: random.Random) -> KMatrix:
    A = [[F.z() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            x = F.rand(rng)
            A[i][j] = x
            A[j][i] = x
    return A


def kidentity(F: GFq3, n: int) -> KMatrix:
    return [[F.o() if i == j else F.z() for j in range(n)] for i in range(n)]


def compact_secret_matrix(F: GFq3, G: KMatrix) -> KMatrix:
    V, O = len(G), len(G[0])
    N = V + O
    S = [[F.z() for _ in range(N)] for _ in range(N)]
    for i in range(V):
        S[i][i] = F.o()
        for j in range(O):
            S[i][V + j] = G[i][j]
    for j in range(O):
        S[V + j][V + j] = F.o()
    return S


def compact_secret_inverse(F: GFq3, G: KMatrix) -> KMatrix:
    V, O = len(G), len(G[0])
    N = V + O
    Sinv = kidentity(F, N)
    for i in range(V):
        for j in range(O):
            Sinv[i][V + j] = F.neg(G[i][j])
    return Sinv


def central_forms(F: GFq3, V: int, O: int, m: int, rng: random.Random) -> List[KMatrix]:
    N = V + O
    out: List[KMatrix] = []
    for _ in range(m):
        C = [[F.z() for _ in range(N)] for _ in range(N)]
        A = random_ksymmetric(F, V, rng)
        E = [[F.rand(rng) for _ in range(O)] for _ in range(V)]
        for i in range(V):
            for j in range(V):
                C[i][j] = A[i][j]
            for j in range(O):
                C[i][V + j] = E[i][j]
                C[V + j][i] = E[i][j]
        out.append(C)
    return out


def public_eval_k(F: GFq3, public_forms: Sequence[KMatrix], x: Vector) -> Vector:
    X = base_to_k(x)
    return [F.trace(kquad(F, A, X)) for A in public_forms]


def kforms_to_base_matrices(F: GFq3, forms: Sequence[KMatrix]) -> List[Matrix]:
    """Recover exact base-field symmetric matrices for Tr(x^T P x)."""
    n = 3 * len(forms[0])
    q = F.q
    inv2 = pow(2, q - 2, q)
    basis = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    out: List[Matrix] = []
    for A in forms:
        diag = [public_eval_k(F, [A], e)[0] for e in basis]
        M = [[0] * n for _ in range(n)]
        for i in range(n):
            M[i][i] = diag[i]
        for i in range(n):
            for j in range(i + 1, n):
                e = [(basis[i][k] + basis[j][k]) % q for k in range(n)]
                val = public_eval_k(F, [A], e)[0]
                c = (val - diag[i] - diag[j]) * inv2 % q
                M[i][j] = M[j][i] = c
        out.append(M)
    return out


@dataclass
class ReducedPublicKey:
    q: int
    V: int
    O: int
    public_matrices: List[Matrix]
    hash_prefix: bytes


@dataclass
class ReducedSecretKey:
    F: GFq3
    G: KMatrix
    Sinv: KMatrix
    central_k: List[KMatrix]
    central_base: List[Matrix]


def reduced_keygen(q: int, V: int, O: int, rng: random.Random) -> Tuple[ReducedPublicKey, ReducedSecretKey]:
    F = GFq3(q)
    m = 3 * O
    G = [[F.rand(rng) for _ in range(O)] for _ in range(V)]
    S = compact_secret_matrix(F, G)
    Sinv = compact_secret_inverse(F, G)
    central = central_forms(F, V, O, m, rng)
    public_k = [kmat_mul(F, ktranspose(S), kmat_mul(F, C, S)) for C in central]
    public_base = kforms_to_base_matrices(F, public_k)
    central_base = kforms_to_base_matrices(F, central)
    prefix = bytes(rng.randrange(256) for _ in range(32))
    pk = ReducedPublicKey(q, V, O, public_base, prefix)
    sk = ReducedSecretKey(F, G, Sinv, central, central_base)
    # Exact public/secret consistency on basis vectors and random points.
    n = 3 * (V + O)
    for _ in range(8):
        x = [rng.randrange(q) for _ in range(n)]
        X = base_to_k(x)
        SX = kmat_vec(F, S, X)
        left = h.eval_map(public_base, x, q)
        right = [F.trace(kquad(F, C, SX)) for C in central]
        assert left == right
    return pk, sk


def hash_to_field(prefix: bytes, message: bytes, salt: bytes, q: int, m: int) -> Vector:
    """Reduced SHAKE256 hash-to-field with unbiased byte rejection."""
    # Generate more bytes on demand by lengthening SHAKE output.
    bound = (256 // q) * q
    out: Vector = []
    length = 64
    while len(out) < m:
        stream = hashlib.shake_256(prefix + message + salt).digest(length)
        out = [b % q for b in stream if b < bound][:m]
        length *= 2
    return out


def solve_linear_system(A: Matrix, b: Vector, q: int) -> Optional[Vector]:
    """Return one solution to A x=b over F_q, or None."""
    M = [row[:] + [bb % q] for row, bb in zip(A, b)]
    nr = len(M)
    nc = len(A[0]) if A else 0
    pivots: List[int] = []
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if M[i][c] % q), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c] % q, q - 2, q)
        M[r] = [(x * inv) % q for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c] % q:
                z = M[i][c] % q
                M[i] = [(a - z * bb) % q for a, bb in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    for row in M:
        if all(row[c] % q == 0 for c in range(nc)) and row[-1] % q:
            return None
    x = [0] * nc
    for i, c in enumerate(pivots):
        x[c] = M[i][-1] % q
    return x


def honest_sign(pk: ReducedPublicKey, sk: ReducedSecretKey, message: bytes, rng: random.Random) -> Tuple[bytes, Vector]:
    q, V, O = pk.q, pk.V, pk.O
    v, m = 3 * V, 3 * O
    n = v + m
    for _ in range(10000):
        salt = bytes(rng.randrange(256) for _ in range(16))
        target = hash_to_field(pk.hash_prefix, message, salt, q, m)
        vinegar = [rng.randrange(q) for _ in range(v)]
        # Central map is affine-linear in oil variables.  Recover the system by
        # evaluation, avoiding any basis convention ambiguity.
        zero_oil = vinegar + [0] * m
        const = h.eval_map(sk.central_base, zero_oil, q)
        A = [[0] * m for _ in range(m)]
        for j in range(m):
            x = zero_oil[:]
            x[v + j] = 1
            col = h.eval_map(sk.central_base, x, q)
            for i in range(m):
                A[i][j] = (col[i] - const[i]) % q
        rhs = [(target[i] - const[i]) % q for i in range(m)]
        oil = solve_linear_system(A, rhs, q)
        if oil is None:
            continue
        central_x = vinegar + oil
        public_x = k_to_base(kmat_vec(sk.F, sk.Sinv, base_to_k(central_x)))
        if verify(pk, message, (salt, public_x)):
            return salt, public_x
    raise RuntimeError("honest signing failed after 10000 trials")


def verify(pk: ReducedPublicKey, message: bytes, signature: Tuple[bytes, Vector]) -> bool:
    salt, x = signature
    m = 3 * pk.O
    if len(x) != 3 * (pk.V + pk.O):
        return False
    target = hash_to_field(pk.hash_prefix, message, salt, pk.q, m)
    return h.eval_map(pk.public_matrices, x, pk.q) == target


# ---------------------------------------------------------------------------
# Affine quadrics and coefficientwise homotopy lifting over F_q
# ---------------------------------------------------------------------------


@dataclass
class Quad:
    # f(x)=c+b^T x+sum_{i,j} A[i][j]x_i x_j.
    A: Matrix
    b: Vector
    c: int


def fp_series_add(a: Sequence[int], b: Sequence[int], p: int, k: int) -> Vector:
    return [((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(k)]


def fp_series_scale(a: Sequence[int], c: int, p: int, k: int) -> Vector:
    return [((a[i] if i < len(a) else 0) * c) % p for i in range(k)]


def fp_series_mul(a: Sequence[int], b: Sequence[int], p: int, k: int, trace: "AttackTrace") -> Vector:
    out = [0] * k
    for i, x in enumerate(a[:k]):
        if not x:
            continue
        for j, y in enumerate(b[: k - i]):
            if y:
                out[i + j] = (out[i + j] + x * y) % p
                trace.base_series_coeff_mults += 1
    trace.base_series_products += 1
    return out


def quad_eval_series(qd: Quad, X: Sequence[Sequence[int]], p: int, k: int, trace: "AttackTrace") -> Vector:
    out = [qd.c % p] + [0] * (k - 1)
    a = len(X)
    for i in range(a):
        out = fp_series_add(out, fp_series_scale(X[i], qd.b[i], p, k), p, k)
        for j in range(a):
            if qd.A[i][j] % p:
                prod = fp_series_mul(X[i], X[j], p, k, trace)
                out = fp_series_add(out, fp_series_scale(prod, qd.A[i][j], p, k), p, k)
    return out


def quad_jacobian_at(qd: Quad, x: Vector, p: int) -> Vector:
    a = len(x)
    row: Vector = []
    for j in range(a):
        s = qd.b[j] % p
        for r in range(a):
            s = (s + (qd.A[j][r] + qd.A[r][j]) * x[r]) % p
        row.append(s)
    return row


def fp_mat_inv(A: Matrix, p: int) -> Matrix:
    n = len(A)
    M = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c] % p), None)
        if piv is None:
            raise ZeroDivisionError("singular matrix")
        M[c], M[piv] = M[piv], M[c]
        inv = pow(M[c][c] % p, p - 2, p)
        M[c] = [(x * inv) % p for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] % p:
                z = M[i][c] % p
                M[i] = [(a - z * b) % p for a, b in zip(M[i], M[c])]
    return [row[n:] for row in M]


def affine_linear_product_quad(alpha: int, beta: int, a: int, p: int) -> Quad:
    la = [pow(alpha, j, p) for j in range(a)]
    lb = [pow(beta, j, p) for j in range(a)]
    ca, cb = pow(alpha, a, p), pow(beta, a, p)
    A = [[la[i] * lb[j] % p for j in range(a)] for i in range(a)]
    b = [(ca * lb[i] + cb * la[i]) % p for i in range(a)]
    return Quad(A, b, ca * cb % p)


def poly_from_roots_fp(roots: Sequence[int], p: int) -> Vector:
    out = [1]
    for r in roots:
        nxt = [0] * (len(out) + 1)
        for i, c in enumerate(out):
            nxt[i] = (nxt[i] - r * c) % p
            nxt[i + 1] = (nxt[i + 1] + c) % p
        out = nxt
    return out


def start_root(selected_nodes: Sequence[int], p: int) -> Vector:
    coeff = poly_from_roots_fp(selected_nodes, p)
    return [coeff[j] % p for j in range(len(selected_nodes))]


def lift_branch(g: Sequence[Quad], f: Sequence[Quad], x0: Vector, p: int, precision: int, trace: "AttackTrace") -> List[Vector]:
    a = len(x0)
    X = [[x0[i]] + [0] * (precision - 1) for i in range(a)]
    J0 = [quad_jacobian_at(gi, x0, p) for gi in g]
    J0inv = fp_mat_inv(J0, p)
    for r in range(1, precision):
        # H=(1-t)g+t f = g + t(f-g).
        residual: Vector = []
        for gi, fi in zip(g, f):
            gv = quad_eval_series(gi, X, p, r + 1, trace)
            fv = quad_eval_series(fi, X, p, r + 1, trace)
            coeff = (gv[r] + (fv[r - 1] - gv[r - 1])) % p
            residual.append(coeff)
        delta = [(-sum(J0inv[i][j] * residual[j] for j in range(a))) % p for i in range(a)]
        for i in range(a):
            X[i][r] = delta[i]
    trace.branch_lifts += 1
    trace.branch_coefficients += a * precision
    return X


def verify_lifted_branch(g: Sequence[Quad], f: Sequence[Quad], X: Sequence[Vector], p: int, precision: int) -> None:
    """Exact regression check H(t,X(t))=0 mod t^precision, off-ledger."""
    dummy = AttackTrace()
    for gi, fi in zip(g, f):
        gv = quad_eval_series(gi, X, p, precision, dummy)
        fv = quad_eval_series(fi, X, p, precision, dummy)
        for k in range(precision):
            coeff = gv[k]
            if k >= 1:
                coeff = (coeff + fv[k - 1] - gv[k - 1]) % p
            if coeff % p:
                raise ArithmeticError(f"lifted branch residual is nonzero at t^{k}")


def quad_from_chart(R: Matrix, z0: Vector, B: Matrix, p: int) -> Quad:
    # z=z0+B x, where columns of B span the chart kernel.
    Bt = h.transpose(B)
    A = h.matmul(Bt, h.matmul(R, B, p), p)
    Rz0 = h.matvec(R, z0, p)
    b = [(2 * x) % p for x in h.matvec(Bt, Rz0, p)]
    c = h.quad(R, z0, p)
    return Quad(A, b, c)


def random_chart(r: int, p: int, rng: random.Random) -> Tuple[Vector, Vector, Matrix]:
    while True:
        ell = [rng.randrange(p) for _ in range(r)]
        if any(ell):
            break
    pivot = next(i for i, x in enumerate(ell) if x)
    inv = pow(ell[pivot], p - 2, p)
    z0 = [0] * r
    z0[pivot] = inv
    cols: List[Vector] = []
    for j in range(r):
        if j == pivot:
            continue
        col = [0] * r
        col[j] = 1
        col[pivot] = (-ell[j] * inv) % p
        cols.append(col)
    B = h.columns_to_matrix(cols, r)
    assert sum(ell[i] * z0[i] for i in range(r)) % p == 1
    assert all(sum(ell[i] * B[i][j] for i in range(r)) % p == 0 for j in range(r - 1))
    return ell, z0, B


# ---------------------------------------------------------------------------
# K' = F_{q^2}, extension-valued series, recombination, reconstruction
# ---------------------------------------------------------------------------


FE = Tuple[int, int]
Series = List[FE]
UPolySeries = List[Series]  # coefficient list in U, low to high
PolyFE = List[FE]  # polynomial in t, low to high


@dataclass
class ExtCounter:
    add: int = 0
    mul: int = 0
    inv: int = 0
    pow: int = 0


class GFp2:
    """F_p[w]/(w^2-nu) with nu a quadratic non-residue."""

    def __init__(self, p: int, counter: Optional[ExtCounter] = None):
        self.p = p
        self.counter = counter or ExtCounter()
        self.nu = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)

    def z(self) -> FE:
        return (0, 0)

    def o(self) -> FE:
        return (1, 0)

    def from_fp(self, x: int) -> FE:
        return (x % self.p, 0)

    def is_zero(self, x: FE) -> bool:
        return x[0] % self.p == 0 and x[1] % self.p == 0

    def add(self, x: FE, y: FE) -> FE:
        self.counter.add += 1
        p = self.p
        return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)

    def neg(self, x: FE) -> FE:
        p = self.p
        return ((-x[0]) % p, (-x[1]) % p)

    def sub(self, x: FE, y: FE) -> FE:
        return self.add(x, self.neg(y))

    def mul(self, x: FE, y: FE) -> FE:
        self.counter.mul += 1
        p = self.p
        a = (x[0] * y[0] + self.nu * x[1] * y[1]) % p
        b = (x[0] * y[1] + x[1] * y[0]) % p
        return (a, b)

    def inv(self, x: FE) -> FE:
        self.counter.inv += 1
        if self.is_zero(x):
            raise ZeroDivisionError
        p = self.p
        den = (x[0] * x[0] - self.nu * x[1] * x[1]) % p
        iden = pow(den, p - 2, p)
        return (x[0] * iden % p, -x[1] * iden % p)

    def div(self, x: FE, y: FE) -> FE:
        return self.mul(x, self.inv(y))

    def pow(self, x: FE, e: int) -> FE:
        self.counter.pow += 1
        out = self.o()
        base = x
        while e:
            if e & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            e >>= 1
        return out

    def frob(self, x: FE) -> FE:
        # For a nonsquare nu, w^p=-w.
        return (x[0] % self.p, (-x[1]) % self.p)

    def rand(self, rng: random.Random) -> FE:
        return (rng.randrange(self.p), rng.randrange(self.p))

    def elements(self) -> Iterable[FE]:
        for a in range(self.p):
            for b in range(self.p):
                yield (a, b)


def s_zero(F: GFp2, P: int) -> Series:
    return [F.z() for _ in range(P)]


def s_const(F: GFp2, x: FE, P: int) -> Series:
    out = s_zero(F, P)
    out[0] = x
    return out


def s_add(F: GFp2, a: Series, b: Series) -> Series:
    return [F.add(x, y) for x, y in zip(a, b)]


def s_neg(F: GFp2, a: Series) -> Series:
    return [F.neg(x) for x in a]


def s_sub(F: GFp2, a: Series, b: Series) -> Series:
    return [F.sub(x, y) for x, y in zip(a, b)]


def s_scale(F: GFp2, a: Series, c: FE) -> Series:
    return [F.mul(x, c) for x in a]


def s_mul(F: GFp2, a: Series, b: Series, trace: "AttackTrace") -> Series:
    P = len(a)
    out = s_zero(F, P)
    for i, x in enumerate(a):
        if F.is_zero(x):
            continue
        for j in range(P - i):
            y = b[j]
            if not F.is_zero(y):
                out[i + j] = F.add(out[i + j], F.mul(x, y))
                trace.extension_series_coeff_mults += 1
    trace.extension_series_products += 1
    return out


def s_inv(F: GFp2, a: Series, trace: "AttackTrace") -> Series:
    if F.is_zero(a[0]):
        raise ZeroDivisionError("series has zero constant coefficient")
    P = len(a)
    out = s_zero(F, P)
    out[0] = F.inv(a[0])
    for n in range(1, P):
        acc = F.z()
        for i in range(1, n + 1):
            acc = F.add(acc, F.mul(a[i], out[n - i]))
            trace.extension_series_coeff_mults += 1
        out[n] = F.neg(F.mul(out[0], acc))
    trace.extension_series_inversions += 1
    return out


def s_div(F: GFp2, a: Series, b: Series, trace: "AttackTrace") -> Series:
    return s_mul(F, a, s_inv(F, b, trace), trace)


def upoly_zero(F: GFp2, degree: int, P: int) -> UPolySeries:
    return [s_zero(F, P) for _ in range(degree + 1)]


def upoly_add(F: GFp2, A: UPolySeries, B: UPolySeries) -> UPolySeries:
    n = max(len(A), len(B))
    P = len(A[0] if A else B[0])
    out = upoly_zero(F, n - 1, P)
    for i in range(n):
        if i < len(A):
            out[i] = s_add(F, out[i], A[i])
        if i < len(B):
            out[i] = s_add(F, out[i], B[i])
    return out


def upoly_mul(F: GFp2, A: UPolySeries, B: UPolySeries, trace: "AttackTrace") -> UPolySeries:
    P = len(A[0])
    out = upoly_zero(F, len(A) + len(B) - 2, P)
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            out[i + j] = s_add(F, out[i + j], s_mul(F, a, b, trace))
    trace.upoly_products += 1
    return out


def product_tree_recombine(
    F: GFp2, u: Sequence[Series], scalar_families: Sequence[Sequence[Series]], trace: "AttackTrace"
) -> Tuple[UPolySeries, List[UPolySeries]]:
    """Use q=qLqR and v=vLqR+qLvR at every product-tree node."""
    P = len(u[0])
    nodes: List[Tuple[UPolySeries, List[UPolySeries]]] = []
    for s in range(len(u)):
        qleaf = [s_neg(F, u[s]), s_const(F, F.o(), P)]
        vleaves = [[scalar_families[k][s]] for k in range(len(scalar_families))]
        nodes.append((qleaf, vleaves))
    while len(nodes) > 1:
        nxt: List[Tuple[UPolySeries, List[UPolySeries]]] = []
        for i in range(0, len(nodes), 2):
            if i + 1 == len(nodes):
                nxt.append(nodes[i])
                continue
            qL, vL = nodes[i]
            qR, vR = nodes[i + 1]
            q = upoly_mul(F, qL, qR, trace)
            vv: List[UPolySeries] = []
            for a, b in zip(vL, vR):
                vv.append(upoly_add(F, upoly_mul(F, a, qR, trace), upoly_mul(F, qL, b, trace)))
            nxt.append((q, vv))
            trace.product_tree_nodes += 1
        nodes = nxt
    return nodes[0]


def fe_linear_solve(F: GFp2, A: List[List[FE]], b: List[FE]) -> Optional[List[FE]]:
    M = [row[:] + [bb] for row, bb in zip(A, b)]
    nr = len(M)
    nc = len(A[0]) if A else 0
    pivots: List[int] = []
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if not F.is_zero(M[i][c])), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = F.inv(M[r][c])
        M[r] = [F.mul(x, inv) for x in M[r]]
        for i in range(nr):
            if i != r and not F.is_zero(M[i][c]):
                z = M[i][c]
                M[i] = [F.sub(x, F.mul(z, y)) for x, y in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    for row in M:
        if all(F.is_zero(row[c]) for c in range(nc)) and not F.is_zero(row[-1]):
            return None
    x = [F.z() for _ in range(nc)]
    for i, c in enumerate(pivots):
        x[c] = M[i][-1]
    return x


def trim_poly(F: GFp2, a: PolyFE) -> PolyFE:
    out = a[:]
    while len(out) > 1 and F.is_zero(out[-1]):
        out.pop()
    return out


def rational_reconstruct(F: GFp2, series: Series, max_num: int, max_den: int, trace: "AttackTrace") -> Tuple[PolyFE, PolyFE]:
    """Find N,D with D(0)=1 and D*S=N mod t^P."""
    P = len(series)
    unknowns = max_den + (max_num + 1)  # d1..dD, n0..nN
    A: List[List[FE]] = []
    rhs: List[FE] = []
    for k in range(P):
        row = [F.z() for _ in range(unknowns)]
        for j in range(1, max_den + 1):
            if k - j >= 0:
                row[j - 1] = series[k - j]
        if k <= max_num:
            row[max_den + k] = F.neg(F.o())
        A.append(row)
        rhs.append(F.neg(series[k]))
    sol = fe_linear_solve(F, A, rhs)
    if sol is None:
        raise ArithmeticError("rational reconstruction inconsistent")
    D = [F.o()] + sol[:max_den]
    N = sol[max_den:]
    # Verify every Taylor coefficient.
    for k in range(P):
        lhs = F.z()
        for j in range(min(k, max_den) + 1):
            lhs = F.add(lhs, F.mul(D[j], series[k - j]))
        rhs_k = N[k] if k <= max_num else F.z()
        if lhs != rhs_k:
            raise ArithmeticError("rational reconstruction verification failed")
    trace.rational_reconstructions += 1
    return trim_poly(F, N), trim_poly(F, D)


def poly_shift_lead_at_one(F: GFp2, p: PolyFE) -> Tuple[int, FE]:
    """Return ord_{t=1} p and the first nonzero coefficient of p(1+s)."""
    if all(F.is_zero(x) for x in p):
        return 10**9, F.z()
    for k in range(len(p)):
        acc = F.z()
        for i in range(k, len(p)):
            coeff = math.comb(i, k) % F.p
            if coeff:
                acc = F.add(acc, F.mul(p[i], F.from_fp(coeff)))
        if not F.is_zero(acc):
            return k, acc
    raise AssertionError("nonzero polynomial had no shifted coefficient")


def rational_valuation_lead(F: GFp2, rat: Tuple[PolyFE, PolyFE]) -> Tuple[int, FE]:
    N, D = rat
    on, ln = poly_shift_lead_at_one(F, N)
    od, ld = poly_shift_lead_at_one(F, D)
    if on >= 10**9:
        return 10**9, F.z()
    return on - od, F.div(ln, ld)


def specialize_upoly(F: GFp2, rats: Sequence[Tuple[PolyFE, PolyFE]], common_val: int) -> List[FE]:
    out: List[FE] = []
    for rat in rats:
        val, lead = rational_valuation_lead(F, rat)
        if val < common_val:
            raise ArithmeticError(f"numerator valuation {val} below q valuation {common_val}")
        out.append(lead if val == common_val else F.z())
    return trim_poly(F, out)


def poly_eval_fe(F: GFp2, p: Sequence[FE], x: FE) -> FE:
    out = F.z()
    for c in reversed(p):
        out = F.add(F.mul(out, x), c)
    return out


def poly_derivative_fe(F: GFp2, p: Sequence[FE]) -> List[FE]:
    if len(p) <= 1:
        return [F.z()]
    return [F.mul(F.from_fp(i), p[i]) for i in range(1, len(p))]


def base_series_to_ext(F: GFp2, a: Sequence[int]) -> Series:
    return [F.from_fp(x) for x in a]


def ext_dot_series(F: GFp2, coeffs: Sequence[FE], X: Sequence[Series]) -> Series:
    P = len(X[0])
    out = s_zero(F, P)
    for c, x in zip(coeffs, X):
        out = s_add(F, out, s_scale(F, x, c))
    return out


def eval_quad_ext_series(F: GFp2, qd: Quad, X: Sequence[Series], trace: "AttackTrace") -> Series:
    P = len(X[0])
    out = s_const(F, F.from_fp(qd.c), P)
    a = len(X)
    for i in range(a):
        if qd.b[i] % F.p:
            out = s_add(F, out, s_scale(F, X[i], F.from_fp(qd.b[i])))
        for j in range(a):
            if qd.A[i][j] % F.p:
                out = s_add(
                    F,
                    out,
                    s_scale(F, s_mul(F, X[i], X[j], trace), F.from_fp(qd.A[i][j])),
                )
    return out


def quad_eval_fe(F: GFp2, qd: Quad, x: Sequence[FE]) -> FE:
    out = F.from_fp(qd.c)
    for i in range(len(x)):
        if qd.b[i] % F.p:
            out = F.add(out, F.mul(F.from_fp(qd.b[i]), x[i]))
        for j in range(len(x)):
            if qd.A[i][j] % F.p:
                out = F.add(out, F.mul(F.from_fp(qd.A[i][j]), F.mul(x[i], x[j])))
    return out


def exhaustive_core_crosscheck(
    p: int,
    f: Sequence[Quad],
    lam: Sequence[FE],
    qbar: Sequence[FE],
    accepted_tau: FE,
    accepted_coords: Vector,
) -> Dict[str, object]:
    """Independent exhaustive check for the tiny a=2 demonstration core."""
    if len(f) > 2:
        return {"skipped": True, "reason": "core dimension exceeds two"}
    E = GFp2(p, ExtCounter())
    base_roots: List[Vector] = []
    for x in itertools.product(range(p), repeat=len(f)):
        xx = list(x)
        if all(quad_eval_fp(qd, xx, p) == 0 for qd in f):
            base_roots.append(xx)
    ext_roots: List[List[FE]] = []
    projections: Dict[FE, int] = {}
    elems = list(E.elements())
    for x in itertools.product(elems, repeat=len(f)):
        if all(E.is_zero(quad_eval_fe(E, qd, x)) for qd in f):
            xx = list(x)
            ext_roots.append(xx)
            tau = E.z()
            for c, xi in zip(lam, xx):
                tau = E.add(tau, E.mul(c, xi))
            if not E.is_zero(poly_eval_fe(E, qbar, tau)):
                raise ArithmeticError("enumerated extension root is absent from target eliminant")
            projections[tau] = projections.get(tau, 0) + 1
    qprime = poly_derivative_fe(E, qbar)
    simple_roots = [x for x in elems if E.is_zero(poly_eval_fe(E, qbar, x)) and not E.is_zero(poly_eval_fe(E, qprime, x))]
    for tau in simple_roots:
        if projections.get(tau, 0) != 1:
            raise ArithmeticError("simple eliminant root does not have exactly one enumerated F_{q^2} fiber point")
    if accepted_coords not in base_roots:
        raise ArithmeticError("accepted coordinates are absent from exhaustive base-field roots")
    if projections.get(accepted_tau, 0) != 1:
        raise ArithmeticError("accepted separator value does not isolate one extension-field point")
    return {
        "skipped": False,
        "base_root_count": len(base_roots),
        "extension_root_count": len(ext_roots),
        "simple_eliminant_roots_in_Fq2": len(simple_roots),
        "accepted_fiber_size": projections.get(accepted_tau, 0),
    }


# ---------------------------------------------------------------------------
# End-to-end local-separator attack
# ---------------------------------------------------------------------------


@dataclass
class AttackTrace:
    branch_lifts: int = 0
    branch_coefficients: int = 0
    base_series_products: int = 0
    base_series_coeff_mults: int = 0
    extension_series_products: int = 0
    extension_series_coeff_mults: int = 0
    extension_series_inversions: int = 0
    upoly_products: int = 0
    product_tree_nodes: int = 0
    rational_reconstructions: int = 0
    coordinate_reruns: int = 0
    target_roots_examined: int = 0
    simple_target_roots: int = 0
    frobenius_passes: int = 0
    filter_passes: int = 0
    exact_candidates_verified: int = 0
    branch_residual_checks: int = 0
    exhaustive_base_roots: int = 0
    exhaustive_extension_roots: int = 0
    simple_fibers_crosschecked: int = 0
    stage_seconds: Dict[str, float] = field(default_factory=dict)
    ext_ops: Dict[str, int] = field(default_factory=dict)


@dataclass
class AttackAttempt:
    success: bool
    failure: str
    seed_points: int
    cw_points: int
    chart_tries: int
    separator_tries: int
    target_eliminant_degree: int
    trace: AttackTrace
    certificate: Optional[Dict[str, object]] = None


@dataclass
class ForgeryResult:
    success: bool
    salt_hex: str
    signature: Vector
    target: Vector
    attempts: int
    honest_signature_verified: bool
    forgery_verified: bool
    parameters: Dict[str, int]
    attempt_records: List[Dict[str, object]]
    public_instance: Dict[str, object]
    successful_certificate: Optional[Dict[str, object]]


def exact_target_transform(y: Vector, p: int) -> Matrix:
    T = h.target_annihilating_transform(y, p)
    Ty = h.matvec(T, y, p)
    if Ty[-1] == 0:
        raise AssertionError
    inv = pow(Ty[-1], p - 2, p)
    T[-1] = [(inv * x) % p for x in T[-1]]
    assert h.matvec(T, y, p) == [0] * (len(y) - 1) + [1]
    return T


def seed_anchor(forms: Sequence[Matrix], p: int, rng: random.Random) -> Tuple[Matrix, Vector, int]:
    n = len(forms[0])
    ambient = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    S0 = h.random_subspace_from_complement(ambient, 7, p, rng)
    restricted = [h.restrict_form(A, S0, p) for A in forms]
    tested = 0
    for z in h.projective_points(7, p):
        tested += 1
        if all(h.quad(A, z, p) == 0 for A in restricted):
            w = h.matvec(S0, z, p)
            return S0, w, tested
    raise AssertionError("Chevalley--Warning anchor not found")


def construct_common_isotropic_anchored(
    forms: Sequence[Matrix], anchor: Vector, r: int, p: int, rng: random.Random
) -> Tuple[Matrix, int]:
    s = len(forms)
    n = len(anchor)
    W_cols = [anchor[:]]
    tested = 0
    for _j in range(1, r):
        constraints: Matrix = []
        for A in forms:
            for w in W_cols:
                constraints.append(h.vecmat(w, A, p))
        Lbasis = h.nullspace(constraints, n, p)
        extended = h.independent_extend(W_cols, Lbasis, p)
        comp = extended[len(W_cols) :]
        d = 2 * s + 1
        if len(comp) < d:
            raise ArithmeticError("Chevalley--Warning quotient dimension too small")
        Smat = h.random_subspace_from_complement(comp, d, p, rng)
        restricted = [h.restrict_form(A, Smat, p) for A in forms]
        root = None
        for z in h.projective_points(d, p):
            tested += 1
            if all(h.quad(A, z, p) == 0 for A in restricted):
                root = z
                break
        if root is None:
            raise AssertionError("Chevalley--Warning extension root not found")
        W_cols.append(h.matvec(Smat, root, p))
        W = h.columns_to_matrix(W_cols, n)
        assert all(all(x == 0 for row in h.restrict_form(A, W, p) for x in row) for A in forms)
    return h.columns_to_matrix(W_cols, n), tested


def find_fp_square_root(a: int, p: int) -> Optional[int]:
    return next((x for x in range(p) if x * x % p == a % p), None)


def local_separator_solve(
    f: Sequence[Quad],
    H0: Quad,
    p: int,
    rng: random.Random,
    max_separator_tries: int = 50,
    forced_lambda: Optional[Sequence[FE]] = None,
    forced_c: Optional[FE] = None,
) -> Tuple[Optional[Vector], str, int, int, AttackTrace, Optional[Dict[str, object]]]:
    a = len(f)
    C = 1 << a
    C0 = (1 << (a - 1)) * (a + 2)
    precision = 4 * C0 + 1
    trace = AttackTrace()
    t0 = time.perf_counter()
    nodes = list(range(1, 2 * a + 1))
    if 2 * a >= p:
        return None, "field too small for distinct Vandermonde nodes", 0, 0, trace, None
    g = [affine_linear_product_quad(nodes[2 * i], nodes[2 * i + 1], a, p) for i in range(a)]
    starts = [start_root([nodes[2 * i + bit] for i, bit in enumerate(choice)], p) for choice in itertools.product((0, 1), repeat=a)]
    branches = [lift_branch(g, f, x0, p, precision, trace) for x0 in starts]
    for branch in branches:
        verify_lifted_branch(g, f, branch, p, precision)
        trace.branch_residual_checks += 1
    trace.stage_seconds["initial_branch_lifts"] = time.perf_counter() - t0

    counter = ExtCounter()
    E = GFp2(p, counter)
    ext_branches = [[base_series_to_ext(E, coord) for coord in branch] for branch in branches]

    sep_budget = 1 if forced_lambda is not None else max_separator_tries
    for sep_try in range(1, sep_budget + 1):
        t1 = time.perf_counter()
        lam = list(forced_lambda) if forced_lambda is not None else [E.rand(rng) for _ in range(a)]
        u = [ext_dot_series(E, lam, X) for X in ext_branches]
        start_values = {us[0] for us in u}
        allowed_c = [x for x in E.elements() if E.neg(x) not in start_values]
        if not allowed_c:
            continue
        if forced_c is not None:
            c = forced_c
            if E.neg(c) in start_values:
                continue
        else:
            c = allowed_c[rng.randrange(len(allowed_c))]
        lam_qinv = [E.frob(x) for x in lam]  # d=2: q^(d-1)=q.
        b1 = [ext_dot_series(E, lam_qinv, X) for X in ext_branches]
        psi: List[Series] = []
        try:
            for X, us in zip(ext_branches, u):
                num = eval_quad_ext_series(E, H0, X, trace)
                den = s_add(E, s_const(E, c, precision), us)
                psi.append(s_div(E, num, den, trace))
        except ZeroDivisionError:
            continue
        qser, (vb1ser, vpsiser) = product_tree_recombine(E, u, [b1, psi], trace)
        trace.stage_seconds["scalar_recombination"] = trace.stage_seconds.get("scalar_recombination", 0.0) + time.perf_counter() - t1

        t2 = time.perf_counter()
        try:
            qrats = [rational_reconstruct(E, s, C0, C0, trace) for s in qser]
            b1rats = [rational_reconstruct(E, s, C0, C0, trace) for s in vb1ser]
            psirats = [rational_reconstruct(E, s, 2 * C0, 2 * C0, trace) for s in vpsiser]
        except ArithmeticError:
            continue
        qvals = [rational_valuation_lead(E, r)[0] for r in qrats]
        qval = min(qvals)
        try:
            qbar = specialize_upoly(E, qrats, qval)
            b1bar = specialize_upoly(E, b1rats, qval)
            psibar = specialize_upoly(E, psirats, qval)
        except ArithmeticError:
            continue
        trace.stage_seconds["rational_reconstruction_specialization"] = trace.stage_seconds.get(
            "rational_reconstruction_specialization", 0.0
        ) + time.perf_counter() - t2

        qprime = poly_derivative_fe(E, qbar)
        target_roots = [x for x in E.elements() if E.is_zero(poly_eval_fe(E, qbar, x))]
        trace.target_roots_examined += len(target_roots)
        for tau in target_roots:
            deriv = poly_eval_fe(E, qprime, tau)
            if E.is_zero(deriv):
                continue
            trace.simple_target_roots += 1
            b1val = E.div(poly_eval_fe(E, b1bar, tau), deriv)
            if E.pow(b1val, p) != tau:
                continue
            trace.frobenius_passes += 1
            psival = E.div(poly_eval_fe(E, psibar, tau), deriv)
            h0val = E.mul(psival, E.add(c, tau))
            if h0val[1] != 0 or h0val[0] == 0 or find_fp_square_root(h0val[0], p) is None:
                continue
            trace.filter_passes += 1

            # Exact one-time coordinate rerun.
            t3 = time.perf_counter()
            rerun_branches = [lift_branch(g, f, x0, p, precision, trace) for x0 in starts]
            trace.coordinate_reruns += 1
            ext_rerun = [[base_series_to_ext(E, coord) for coord in branch] for branch in rerun_branches]
            u2 = [ext_dot_series(E, lam, X) for X in ext_rerun]
            coords: Vector = []
            try:
                for j in range(a):
                    family = [X[j] for X in ext_rerun]
                    q2, (vj,) = product_tree_recombine(E, u2, [family], trace)
                    # q2 is recomputed intentionally as part of the literal rerun.
                    q2rats = [rational_reconstruct(E, s, C0, C0, trace) for s in q2]
                    q2val = min(rational_valuation_lead(E, r)[0] for r in q2rats)
                    q2bar = specialize_upoly(E, q2rats, q2val)
                    if q2bar != qbar:
                        raise ArithmeticError("coordinate rerun eliminant mismatch")
                    vjrats = [rational_reconstruct(E, s, C0, C0, trace) for s in vj]
                    vjbar = specialize_upoly(E, vjrats, qval)
                    xj = E.div(poly_eval_fe(E, vjbar, tau), deriv)
                    if xj[1] != 0:
                        raise ArithmeticError("Frobenius-passing point did not descend coordinatewise")
                    coords.append(xj[0])
            except ArithmeticError:
                continue
            trace.stage_seconds["coordinate_rerun"] = trace.stage_seconds.get("coordinate_rerun", 0.0) + time.perf_counter() - t3
            if any(quad_eval_fp(qd, coords, p) != 0 for qd in f):
                continue
            if quad_eval_fp(H0, coords, p) != h0val[0]:
                continue
            trace.exact_candidates_verified += 1
            cross = exhaustive_core_crosscheck(p, f, lam, qbar, tau, coords)
            if not cross.get("skipped"):
                trace.exhaustive_base_roots = int(cross["base_root_count"])
                trace.exhaustive_extension_roots = int(cross["extension_root_count"])
                trace.simple_fibers_crosschecked = int(cross["simple_eliminant_roots_in_Fq2"])
            certificate: Dict[str, object] = {
                "a": a,
                "C": C,
                "C0": C0,
                "precision": precision,
                "vandermonde_nodes": nodes,
                "start_roots": starts,
                "lambda": [list(x) for x in lam],
                "c": list(c),
                "q_valuation_at_target": qval,
                "qbar": [list(x) for x in qbar],
                "b1bar": [list(x) for x in b1bar],
                "psibar": [list(x) for x in psibar],
                "tau": list(tau),
                "qprime_at_tau": list(deriv),
                "b1_at_tau": list(b1val),
                "psi_at_tau": list(psival),
                "H0_at_point": list(h0val),
                "affine_coordinates": coords,
                "crosscheck": cross,
            }
            trace.ext_ops = asdict(counter)
            return coords, "success", sep_try, len(qbar) - 1, trace, certificate
    trace.ext_ops = asdict(counter)
    return None, "no separator produced an accepted simple base-field root", sep_budget, 0, trace, None


def quad_eval_fp(qd: Quad, x: Vector, p: int) -> int:
    out = qd.c % p
    for i in range(len(x)):
        out = (out + qd.b[i] * x[i]) % p
        for j in range(len(x)):
            out = (out + qd.A[i][j] * x[i] * x[j]) % p
    return out


def attack_one_target(
    pk: ReducedPublicKey,
    message: bytes,
    salt: bytes,
    rng: random.Random,
    max_outer_attempts: int,
    max_chart_tries: int,
    max_separator_tries: int,
) -> Tuple[Optional[Vector], List[AttackAttempt], Vector]:
    p = pk.q
    forms = pk.public_matrices
    m = len(forms)
    n = len(forms[0])
    target = hash_to_field(pk.hash_prefix, message, salt, p, m)
    if not any(target):
        return [0] * n, [], target
    T = exact_target_transform(target, p)
    mixed = h.output_mix(forms, T, p)
    ty = h.matvec(T, target, p)
    assert ty == [0] * (m - 1) + [1]
    r = m - 3
    records: List[AttackAttempt] = []

    for _outer in range(max_outer_attempts):
        S0, w, seed_points = seed_anchor(mixed[:3], p, rng)
        U, cw_points = construct_common_isotropic_anchored(mixed[:3], w, r, p, rng)
        restricted = [h.restrict_form(A, U, p) for A in mixed]
        for chart_try in range(1, max_chart_tries + 1):
            _ell, z0, B = random_chart(r, p, rng)
            f = [quad_from_chart(A, z0, B, p) for A in restricted[3:-1]]
            H0 = quad_from_chart(restricted[-1], z0, B, p)
            coords, reason, septries, elimdeg, trace, local_cert = local_separator_solve(
                f, H0, p, rng, max_separator_tries=max_separator_tries
            )
            if coords is None:
                records.append(AttackAttempt(False, reason, seed_points, cw_points, chart_try, septries, elimdeg, trace, None))
                continue
            z = [(z0[i] + sum(B[i][j] * coords[j] for j in range(len(coords)))) % p for i in range(r)]
            h0 = h.quad(restricted[-1], z, p)
            root = find_fp_square_root(pow(h0, p - 2, p), p) if h0 else None
            if root is None:
                records.append(AttackAttempt(False, "accepted filter value had no affine rescaling", seed_points, cw_points, chart_try, septries, elimdeg, trace, local_cert))
                continue
            zscaled = [(root * x) % p for x in z]
            x = h.matvec(U, zscaled, p)
            if h.eval_map(mixed, x, p) != ty:
                records.append(AttackAttempt(False, "candidate failed mixed public map", seed_points, cw_points, chart_try, septries, elimdeg, trace, local_cert))
                continue
            if h.eval_map(forms, x, p) != target:
                records.append(AttackAttempt(False, "candidate failed original public map", seed_points, cw_points, chart_try, septries, elimdeg, trace, local_cert))
                continue
            full_cert: Dict[str, object] = {
                "target_transform": T,
                "isotropic_subspace_U": U,
                "chart_z0": z0,
                "chart_kernel_B": B,
                "affine_system": [asdict(qd) for qd in f],
                "scale_form_H0": asdict(H0),
                "local_separator": local_cert,
                "projective_chart_point": z,
                "affine_scale": root,
                "forged_signature_vector": x,
            }
            records.append(AttackAttempt(True, "success", seed_points, cw_points, chart_try, septries, elimdeg, trace, full_cert))
            return x, records, target
    return None, records, target


def run_forgery(
    q: int,
    V: int,
    O: int,
    seed: int,
    message: bytes,
    max_salts: int,
    max_outer_attempts: int,
    max_chart_tries: int,
    max_separator_tries: int,
) -> ForgeryResult:
    rng = random.Random(seed)
    pk, sk = reduced_keygen(q, V, O, rng)
    public_instance: Dict[str, object] = {
        "q": q,
        "V": V,
        "O": O,
        "hash_prefix_hex": pk.hash_prefix.hex(),
        "message_hex": message.hex(),
        "public_matrices": pk.public_matrices,
    }
    honest = honest_sign(pk, sk, message, rng)
    honest_ok = verify(pk, message, honest)
    records: List[Dict[str, object]] = []
    for salt_index in range(1, max_salts + 1):
        salt = hashlib.sha256(f"attack-salt-{seed}-{salt_index}".encode()).digest()[:16]
        x, attempt_records, target = attack_one_target(
            pk,
            message,
            salt,
            rng,
            max_outer_attempts=max_outer_attempts,
            max_chart_tries=max_chart_tries,
            max_separator_tries=max_separator_tries,
        )
        records.extend(asdict(r) for r in attempt_records)
        if x is not None and verify(pk, message, (salt, x)):
            return ForgeryResult(
                True,
                salt.hex(),
                x,
                target,
                salt_index,
                honest_ok,
                True,
                {"q": q, "ell": 3, "V": V, "O": O, "v": 3 * V, "m": 3 * O, "n": 3 * (V + O), "a": 3 * O - 4},
                records,
                public_instance,
                records[-1].get("certificate") if records else None,
            )
    return ForgeryResult(
        False,
        "",
        [],
        [],
        max_salts,
        honest_ok,
        False,
        {"q": q, "ell": 3, "V": V, "O": O, "v": 3 * V, "m": 3 * O, "n": 3 * (V + O), "a": 3 * O - 4},
        records,
        public_instance,
        None,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=13)
    ap.add_argument("--V", type=int, default=3)
    ap.add_argument("--O", type=int, default=2)
    ap.add_argument("--seed", type=int, default=20260731)
    ap.add_argument("--message", type=str, default="local-separator reduced-parameter forgery")
    ap.add_argument("--max-salts", type=int, default=8)
    ap.add_argument("--max-outer-attempts", type=int, default=8)
    ap.add_argument("--max-chart-tries", type=int, default=4)
    ap.add_argument("--max-separator-tries", type=int, default=30)
    ap.add_argument("--json", type=Path, default=Path("qruov_local_separator_forgery_run.json"))
    args = ap.parse_args()
    if 3 * args.O < 5:
        raise SystemExit("need m=3O at least 5")
    a = 3 * args.O - 4
    if 2 * a >= args.q:
        raise SystemExit("need q>2a for the Vandermonde product start")
    n = 3 * (args.V + args.O)
    r = 3 * args.O - 3
    if n < 4 * r + 3:
        raise SystemExit(f"Chevalley--Warning profile fails n>=4r+3: {n} < {4*r+3}")
    t0 = time.perf_counter()
    result = run_forgery(
        args.q,
        args.V,
        args.O,
        args.seed,
        args.message.encode(),
        args.max_salts,
        args.max_outer_attempts,
        args.max_chart_tries,
        args.max_separator_tries,
    )
    elapsed = time.perf_counter() - t0
    payload = asdict(result)
    payload["elapsed_seconds"] = elapsed
    args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("Reduced QR-UOV local-separator forgery")
    print("parameters:", result.parameters)
    print(f"honest signature verified: {result.honest_signature_verified}")
    print(f"forgery success: {result.success}; verified: {result.forgery_verified}; salt attempts: {result.attempts}")
    print(f"attack records: {len(result.attempt_records)}; elapsed: {elapsed:.3f}s")
    if result.success:
        tr = result.attempt_records[-1]["trace"]
        print("successful trace:", json.dumps(tr, sort_keys=True))
        print("signature salt:", result.salt_hex)
        print("signature vector:", result.signature)
        print("target:", result.target)
    print(f"wrote {args.json}")
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
