#!/usr/bin/env python3
"""Exact toy audit for target-orthogonal Chevalley-Warning slicing.

For a homogeneous quadratic map F_q^n -> F_q^m and a nonzero target y,
we output-mix so that the first three target coordinates vanish. We then
construct an (m-3)-dimensional subspace U on which those three quadrics vanish
identically, using only linear algebra and exhaustive solution of three
quadrics in seven variables. Finally we brute-force the residual square system
on U and verify any returned preimage against the original map.

The script supports random generic MQ maps and hidden-UOV maps. It is intended
as a correctness and success-rate audit on tractable profiles, not a production
benchmark.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from dataclasses import dataclass, asdict
from typing import Iterable, List, Optional, Sequence, Tuple

Vector = List[int]
Matrix = List[List[int]]


def inv_mod(a: int, q: int) -> int:
    a %= q
    if a == 0:
        raise ZeroDivisionError("inverse of zero")
    return pow(a, q - 2, q)


def mat_shape(A: Matrix) -> Tuple[int, int]:
    return (len(A), len(A[0]) if A else 0)


def transpose(A: Matrix) -> Matrix:
    if not A:
        return []
    return [list(row) for row in zip(*A)]


def matmul(A: Matrix, B: Matrix, q: int) -> Matrix:
    if not A or not B:
        return [[] for _ in range(len(A))]
    assert len(A[0]) == len(B)
    Bt = transpose(B)
    return [[sum(x * y for x, y in zip(row, col)) % q for col in Bt] for row in A]


def matvec(A: Matrix, x: Vector, q: int) -> Vector:
    return [sum(a * b for a, b in zip(row, x)) % q for row in A]


def vecmat(x: Vector, A: Matrix, q: int) -> Vector:
    return matvec(transpose(A), x, q)


def dot(x: Vector, y: Vector, q: int) -> int:
    return sum(a * b for a, b in zip(x, y)) % q


def vec_add(x: Vector, y: Vector, q: int) -> Vector:
    return [(a + b) % q for a, b in zip(x, y)]


def vec_scale(a: int, x: Vector, q: int) -> Vector:
    return [(a * v) % q for v in x]


def eye(n: int) -> Matrix:
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def rref(A: Matrix, q: int) -> Tuple[Matrix, List[int]]:
    M = [[v % q for v in row] for row in A]
    if not M:
        return M, []
    rows, cols = len(M), len(M[0])
    pivots: List[int] = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if M[i][c] % q), None)
        if pivot is None:
            continue
        M[r], M[pivot] = M[pivot], M[r]
        z = inv_mod(M[r][c], q)
        M[r] = [(z * v) % q for v in M[r]]
        for i in range(rows):
            if i != r and M[i][c] % q:
                f = M[i][c] % q
                M[i] = [(a - f * b) % q for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return M, pivots


def rank(A: Matrix, q: int) -> int:
    return len(rref(A, q)[1]) if A else 0


def nullspace(A: Matrix, ncols: int, q: int) -> List[Vector]:
    if not A:
        return [row[:] for row in eye(ncols)]
    R, pivots = rref(A, q)
    free = [c for c in range(ncols) if c not in pivots]
    basis: List[Vector] = []
    for f in free:
        x = [0] * ncols
        x[f] = 1
        for i, p in enumerate(pivots):
            x[p] = (-R[i][f]) % q
        basis.append(x)
    return basis


def columns_to_matrix(cols: Sequence[Vector], nrows: int) -> Matrix:
    if not cols:
        return [[] for _ in range(nrows)]
    return [[cols[j][i] for j in range(len(cols))] for i in range(nrows)]


def matrix_columns(A: Matrix) -> List[Vector]:
    return transpose(A) if A and A[0] else []


def independent_extend(seed_cols: Sequence[Vector], candidate_cols: Sequence[Vector], q: int) -> List[Vector]:
    out = [c[:] for c in seed_cols]
    current = rank(columns_to_matrix(out, len(out[0])) if out else [], q) if out else 0
    for c in candidate_cols:
        test = out + [c]
        rr = rank(columns_to_matrix(test, len(c)), q)
        if rr > current:
            out.append(c[:])
            current = rr
    return out


def random_invertible(n: int, q: int, rng: random.Random) -> Matrix:
    while True:
        A = [[rng.randrange(q) for _ in range(n)] for _ in range(n)]
        if rank(A, q) == n:
            return A


def random_symmetric(n: int, q: int, rng: random.Random) -> Matrix:
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = rng.randrange(q)
            A[i][j] = v
            A[j][i] = v
    return A


def quad(A: Matrix, x: Vector, q: int) -> int:
    return dot(x, matvec(A, x, q), q)


def restrict_form(A: Matrix, S: Matrix, q: int) -> Matrix:
    return matmul(transpose(S), matmul(A, S, q), q)


def congruence(A: Matrix, S: Matrix, q: int) -> Matrix:
    return restrict_form(A, S, q)


def output_mix(forms: Sequence[Matrix], T: Matrix, q: int) -> List[Matrix]:
    m = len(forms)
    n = len(forms[0])
    out: List[Matrix] = []
    for row in T:
        A = [[0] * n for _ in range(n)]
        for coeff, F in zip(row, forms):
            if coeff:
                for i in range(n):
                    Ai = A[i]
                    Fi = F[i]
                    for j in range(n):
                        Ai[j] = (Ai[j] + coeff * Fi[j]) % q
        out.append(A)
    assert len(out) == m
    return out


def target_annihilating_transform(y: Vector, q: int) -> Matrix:
    """Return T invertible with T*y = (0,...,0,c), c != 0 for y != 0."""
    m = len(y)
    if all(v % q == 0 for v in y):
        return eye(m)
    perp = nullspace([y], m, q)  # row vectors lambda with lambda dot y = 0
    assert len(perp) == m - 1
    z = next(eye(m)[i] for i, yi in enumerate(y) if yi % q)
    rows = perp + [z]
    assert rank(rows, q) == m
    Ty = matvec(rows, y, q)
    assert all(v == 0 for v in Ty[:-1]) and Ty[-1] != 0
    return rows


def projective_points(dim: int, q: int) -> Iterable[Vector]:
    """Canonical representatives of P^{dim-1}(F_q): first nonzero coordinate is 1."""
    for first in range(dim):
        tail_len = dim - first - 1
        for tail in itertools.product(range(q), repeat=tail_len):
            yield [0] * first + [1] + list(tail)


def random_subspace_from_complement(comp: Sequence[Vector], d: int, q: int, rng: random.Random) -> Matrix:
    assert len(comp) >= d
    ambient = len(comp[0])
    C = columns_to_matrix(comp, ambient)
    while True:
        coeff = [[rng.randrange(q) for _ in range(d)] for _ in range(len(comp))]
        S = matmul(C, coeff, q)
        if rank(S, q) == d:
            return S


@dataclass
class ConstructionStats:
    steps: int = 0
    projective_points_tested: int = 0
    min_quotient_dim: int = 10**9
    max_quotient_dim: int = 0


def construct_common_isotropic(
    forms: Sequence[Matrix], r: int, q: int, rng: random.Random, local_dim: Optional[int] = None
) -> Tuple[Matrix, ConstructionStats]:
    s = len(forms)
    n = len(forms[0])
    d = local_dim if local_dim is not None else 2 * s + 1
    if d <= 2 * s:
        raise ValueError("Chevalley-Warning requires local_dim > sum of degrees = 2s")
    W_cols: List[Vector] = []
    stats = ConstructionStats()
    for j in range(r):
        constraints: Matrix = []
        for A in forms:
            for w in W_cols:
                constraints.append(vecmat(w, A, q))
        L_basis = nullspace(constraints, n, q)
        extended = independent_extend(W_cols, L_basis, q)
        comp = extended[len(W_cols):]
        quotient_dim = len(comp)
        stats.min_quotient_dim = min(stats.min_quotient_dim, quotient_dim)
        stats.max_quotient_dim = max(stats.max_quotient_dim, quotient_dim)
        if quotient_dim < d:
            raise RuntimeError(
                f"dimension bound failed at j={j}: quotient has {quotient_dim}, need {d}"
            )
        S = random_subspace_from_complement(comp, d, q, rng)
        restricted = [restrict_form(A, S, q) for A in forms]
        roots: List[Vector] = []
        tested = 0
        for z in projective_points(d, q):
            tested += 1
            if all(quad(A, z, q) == 0 for A in restricted):
                roots.append(z)
        stats.projective_points_tested += tested
        if not roots:
            raise RuntimeError("Chevalley-Warning audit failed: no nonzero common zero")
        z = roots[rng.randrange(len(roots))]
        v = matvec(S, z, q)
        assert any(v)
        assert all(quad(A, v, q) == 0 for A in forms)
        for A in forms:
            for w in W_cols:
                assert dot(w, matvec(A, v, q), q) == 0
        W_cols.append(v)
        assert rank(columns_to_matrix(W_cols, n), q) == len(W_cols)
        # Exact invariant: every restriction matrix is zero.
        W = columns_to_matrix(W_cols, n)
        assert all(all(x == 0 for row in restrict_form(A, W, q) for x in row) for A in forms)
        stats.steps += 1
    if stats.min_quotient_dim == 10**9:
        stats.min_quotient_dim = n
    return columns_to_matrix(W_cols, n), stats


def brute_force_residual(forms: Sequence[Matrix], target: Vector, U: Matrix, q: int) -> Tuple[Optional[Vector], int]:
    r = len(U[0]) if U and U[0] else 0
    restricted = [restrict_form(A, U, q) for A in forms]
    tested = 0
    for z_tuple in itertools.product(range(q), repeat=r):
        tested += 1
        z = list(z_tuple)
        if all(quad(A, z, q) == t for A, t in zip(restricted, target)):
            return z, tested
    return None, tested


def eval_map(forms: Sequence[Matrix], x: Vector, q: int) -> Vector:
    return [quad(A, x, q) for A in forms]


def generate_generic(n: int, m: int, q: int, rng: random.Random) -> Tuple[List[Matrix], Optional[Matrix]]:
    return [random_symmetric(n, q, rng) for _ in range(m)], None


def generate_uov(n: int, m: int, q: int, rng: random.Random) -> Tuple[List[Matrix], Matrix]:
    if n <= m:
        raise ValueError("UOV requires n > m")
    v, o = n - m, m
    S = random_invertible(n, q, rng)
    forms: List[Matrix] = []
    for _ in range(m):
        C = [[0] * n for _ in range(n)]
        A = random_symmetric(v, q, rng)
        B = [[rng.randrange(q) for _ in range(o)] for _ in range(v)]
        for i in range(v):
            for j in range(v):
                C[i][j] = A[i][j]
            for j in range(o):
                C[i][v + j] = B[i][j]
                C[v + j][i] = B[i][j]
        forms.append(congruence(C, S, q))
    return forms, S


@dataclass
class TrialResult:
    map_kind: str
    q: int
    n: int
    m: int
    r: int
    attempts: int
    success: bool
    residual_points_tested: int
    cw_projective_points_tested: int
    min_quotient_dim: int
    verified: bool
    min_uov_vinegar_projection_rank: Optional[int] = None


def run_target_trial(
    map_kind: str, q: int, n: int, m: int, max_attempts: int, rng: random.Random
) -> TrialResult:
    if map_kind == "generic":
        forms, secret_S = generate_generic(n, m, q, rng)
    elif map_kind == "uov":
        forms, secret_S = generate_uov(n, m, q, rng)
    else:
        raise ValueError(map_kind)
    # Draw a nonzero actual preimage so the original target is guaranteed solvable.
    while True:
        x_star = [rng.randrange(q) for _ in range(n)]
        y = eval_map(forms, x_star, q)
        if any(y):
            break
    T = target_annihilating_transform(y, q)
    mixed = output_mix(forms, T, q)
    ty = matvec(T, y, q)
    assert all(v == 0 for v in ty[:3])
    r = m - 3
    bound = (n - 3) // 4
    if bound < r:
        raise ValueError(f"profile does not support r=m-3: floor((n-3)/4)={bound}, r={r}")

    total_residual = 0
    total_cw = 0
    min_qdim = 10**9
    min_vproj: Optional[int] = None
    for attempt in range(1, max_attempts + 1):
        U, stats = construct_common_isotropic(mixed[:3], r, q, rng)
        if map_kind == "uov":
            assert secret_S is not None
            central_U = matmul(secret_S, U, q)
            vinegar_rows = central_U[: n - m]
            vrank = rank(vinegar_rows, q)
            min_vproj = vrank if min_vproj is None else min(min_vproj, vrank)
        total_cw += stats.projective_points_tested
        min_qdim = min(min_qdim, stats.min_quotient_dim)
        z, tested = brute_force_residual(mixed[3:], ty[3:], U, q)
        total_residual += tested
        if z is None:
            continue
        x = matvec(U, z, q)
        verified = eval_map(forms, x, q) == y
        if not verified:
            raise AssertionError("residual solution failed original-map verification")
        return TrialResult(
            map_kind, q, n, m, r, attempt, True, total_residual, total_cw, min_qdim, True, min_vproj
        )
    return TrialResult(
        map_kind, q, n, m, r, max_attempts, False, total_residual, total_cw, min_qdim, False, min_vproj
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--m", type=int, default=5)
    ap.add_argument("--n", type=int, default=None, help="default 4m-6")
    ap.add_argument("--trials", type=int, default=50)
    ap.add_argument("--max-attempts", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260729)
    ap.add_argument("--kind", choices=["generic", "uov", "both"], default="both")
    ap.add_argument("--json", type=str, default=None)
    args = ap.parse_args()
    n = args.n if args.n is not None else 4 * args.m - 6
    if args.q < 3 or args.q % 2 == 0:
        raise SystemExit("this audit expects an odd prime q")
    rng = random.Random(args.seed)
    kinds = ["generic", "uov"] if args.kind == "both" else [args.kind]
    results: List[TrialResult] = []
    for kind in kinds:
        for _ in range(args.trials):
            results.append(run_target_trial(kind, args.q, n, args.m, args.max_attempts, rng))

    print("Target-orthogonal Chevalley-Warning triplet audit")
    print(f"q={args.q}, n={n}, m={args.m}, r=m-3={args.m-3}, trials/kind={args.trials}")
    print(f"theorem bound floor((n-3)/4)={(n-3)//4}")
    for kind in kinds:
        rr = [x for x in results if x.map_kind == kind]
        succ = sum(x.success for x in rr)
        attempts = [x.attempts for x in rr if x.success]
        first = sum(x.success and x.attempts == 1 for x in rr)
        print(
            f"{kind}: success {succ}/{len(rr)}; first-attempt {first}/{len(rr)}; "
            f"mean attempts among successes={sum(attempts)/len(attempts) if attempts else float('nan'):.3f}; "
            f"all returned candidates verified={all((not x.success) or x.verified for x in rr)}"
        )
        extra = ""
        if kind == "uov":
            vals = [x.min_uov_vinegar_projection_rank for x in rr if x.min_uov_vinegar_projection_rank is not None]
            extra = f"; min UOV vinegar-projection rank={min(vals) if vals else 'n/a'} (full={args.m-3})"
        print(
            f"  min observed extension quotient dimension={min(x.min_quotient_dim for x in rr)}; "
            f"CW projective points tested total={sum(x.cw_projective_points_tested for x in rr)}; "
            f"residual points tested total={sum(x.residual_points_tested for x in rr)}" + extra
        )
    payload = {
        "parameters": vars(args) | {"n_effective": n, "r": args.m - 3, "bound": (n - 3) // 4},
        "results": [asdict(x) for x in results],
    }
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
