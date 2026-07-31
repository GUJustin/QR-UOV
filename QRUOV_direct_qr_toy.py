#!/usr/bin/env python3
"""Exact degree-3 QR-UOV toy for the anchored target-orthogonal reduction.

The toy works over K=F_{q^3}, pulls K-valued UOV quadrics back to F_q by
an explicit trace, hides the oil space with a random K-linear change of
variables, and then performs the public attack entirely over F_q:
  * output-mix so that every target coordinate except one is zero;
  * choose a random seven-dimensional seed subspace and use Chevalley--Warning
    to find a nonzero anchor on the three public quadrics;
  * extend it to an (m-3)-dimensional common totally isotropic subspace
    using the Chevalley-Warning construction;
  * quotient the residual by projective scaling, solve one fewer quadratic,
    and recover the affine scale from the remaining nonzero target equation;
  * verify the returned preimage against the original public map.

This is a correctness regression test, not a production benchmark.
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import random
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

# Reuse audited base-field linear algebra and CW routines.
_HELPER = Path(__file__).resolve().with_name('QRUOV_cw_helper.py')
spec = importlib.util.spec_from_file_location('cwhelper', _HELPER)
h = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = h
assert spec.loader is not None
spec.loader.exec_module(h)

Vector = List[int]
Matrix = List[List[int]]
Elt = Tuple[int, int, int]


class GFq3:
    """F_q[t]/(t^3-t-1), for primes q where this polynomial is irreducible."""

    def __init__(self, q: int):
        self.q = q
        # Verify no root; a cubic is irreducible iff it has no base-field root.
        if any((a**3 - a - 1) % q == 0 for a in range(q)):
            raise ValueError(f't^3-t-1 is reducible over F_{q}')

    def z(self) -> Elt:
        return (0, 0, 0)

    def o(self) -> Elt:
        return (1, 0, 0)

    def add(self, a: Elt, b: Elt) -> Elt:
        q = self.q
        return tuple((a[i] + b[i]) % q for i in range(3))  # type: ignore[return-value]

    def neg(self, a: Elt) -> Elt:
        q = self.q
        return tuple((-x) % q for x in a)  # type: ignore[return-value]

    def sub(self, a: Elt, b: Elt) -> Elt:
        return self.add(a, self.neg(b))

    def mul(self, a: Elt, b: Elt) -> Elt:
        q = self.q
        c = [0] * 5
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i + j] = (c[i + j] + x * y) % q
        # t^3=t+1; t^4=t^2+t. Reduce high-to-low.
        for d in range(4, 2, -1):
            z = c[d] % q
            if z:
                c[d] = 0
                c[d - 2] = (c[d - 2] + z) % q  # t^d -> t^(d-2)+t^(d-3)
                c[d - 3] = (c[d - 3] + z) % q
        return (c[0] % q, c[1] % q, c[2] % q)

    def pow(self, a: Elt, e: int) -> Elt:
        out = self.o()
        base = a
        while e:
            if e & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            e >>= 1
        return out

    def inv(self, a: Elt) -> Elt:
        if a == self.z():
            raise ZeroDivisionError
        return self.pow(a, self.q**3 - 2)

    def rand(self, rng: random.Random) -> Elt:
        return (rng.randrange(self.q), rng.randrange(self.q), rng.randrange(self.q))

    def trace(self, a: Elt) -> int:
        aq = self.pow(a, self.q)
        aq2 = self.pow(aq, self.q)
        return self.add(self.add(a, aq), aq2)[0] % self.q


KMatrix = List[List[Elt]]
KVector = List[Elt]


def kmat_vec(F: GFq3, A: KMatrix, x: KVector) -> KVector:
    out: KVector = []
    for row in A:
        s = F.z()
        for a, b in zip(row, x):
            s = F.add(s, F.mul(a, b))
        out.append(s)
    return out


def kdot(F: GFq3, x: KVector, y: KVector) -> Elt:
    s = F.z()
    for a, b in zip(x, y):
        s = F.add(s, F.mul(a, b))
    return s


def kquad(F: GFq3, A: KMatrix, x: KVector) -> Elt:
    return kdot(F, x, kmat_vec(F, A, x))


def kmat_mul(F: GFq3, A: KMatrix, B: KMatrix) -> KMatrix:
    Bt = [list(col) for col in zip(*B)]
    return [[kdot(F, row, col) for col in Bt] for row in A]


def ktranspose(A: KMatrix) -> KMatrix:
    return [list(col) for col in zip(*A)]


def krank(F: GFq3, A: KMatrix) -> int:
    M = [[x for x in row] for row in A]
    if not M:
        return 0
    nr, nc = len(M), len(M[0])
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if M[i][c] != F.z()), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        z = F.inv(M[r][c])
        M[r] = [F.mul(z, a) for a in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != F.z():
                z = M[i][c]
                M[i] = [F.sub(a, F.mul(z, b)) for a, b in zip(M[i], M[r])]
        r += 1
        if r == nr:
            break
    return r


def random_kinvertible(F: GFq3, n: int, rng: random.Random) -> KMatrix:
    while True:
        A = [[F.rand(rng) for _ in range(n)] for _ in range(n)]
        if krank(F, A) == n:
            return A


def random_ksymmetric(F: GFq3, n: int, rng: random.Random) -> KMatrix:
    A = [[F.z() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            a = F.rand(rng)
            A[i][j] = a
            A[j][i] = a
    return A


def base_to_k(x: Vector) -> KVector:
    assert len(x) % 3 == 0
    return [tuple(x[3*i:3*i+3]) for i in range(len(x)//3)]  # type: ignore[list-item]


def k_to_base(x: KVector) -> Vector:
    return [a for z in x for a in z]


def public_eval(F: GFq3, central_forms: Sequence[KMatrix], S: KMatrix, x: Vector) -> Vector:
    X = base_to_k(x)
    SX = kmat_vec(F, S, X)
    return [F.trace(kquad(F, A, SX)) for A in central_forms]


def public_base_matrices(F: GFq3, central_forms: Sequence[KMatrix], S: KMatrix) -> List[Matrix]:
    n = 3 * len(S)
    q = F.q
    inv2 = pow(2, q - 2, q)
    basis = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    vals = [[public_eval(F, [A], S, e)[0] for e in basis] for A in central_forms]
    out: List[Matrix] = []
    for form_idx, A in enumerate(central_forms):
        M = [[0] * n for _ in range(n)]
        for i in range(n):
            M[i][i] = vals[form_idx][i]
        for i in range(n):
            for j in range(i + 1, n):
                eij = [(basis[i][k] + basis[j][k]) % q for k in range(n)]
                qij = public_eval(F, [A], S, eij)[0]
                mij = (qij - vals[form_idx][i] - vals[form_idx][j]) * inv2 % q
                M[i][j] = M[j][i] = mij
        # Full basis regression.
        for e, v in zip(basis, vals[form_idx]):
            assert h.quad(M, e, q) == v
        out.append(M)
    return out


def generate_qruov(F: GFq3, V: int, O: int, m: int, rng: random.Random):
    N = V + O
    forms: List[KMatrix] = []
    for _ in range(m):
        C = [[F.z() for _ in range(N)] for _ in range(N)]
        A = random_ksymmetric(F, V, rng)
        B = [[F.rand(rng) for _ in range(O)] for _ in range(V)]
        for i in range(V):
            for j in range(V):
                C[i][j] = A[i][j]
            for j in range(O):
                C[i][V + j] = B[i][j]
                C[V + j][i] = B[i][j]
        forms.append(C)
    S = random_kinvertible(F, N, rng)
    mats = public_base_matrices(F, forms, S)
    return forms, S, mats


def seed_anchor(forms: Sequence[Matrix], q: int, rng: random.Random):
    """Choose a random 7-space and find a nonzero common zero of 3 quadrics."""
    if len(forms) != 3:
        raise ValueError('the frozen attack uses three target-zero quadrics')
    n = len(forms[0])
    ambient_basis = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    S0 = h.random_subspace_from_complement(ambient_basis, 7, q, rng)
    restricted = [h.restrict_form(A, S0, q) for A in forms]
    tested = 0
    root = None
    for z in h.projective_points(7, q):
        tested += 1
        if all(h.quad(A, z, q) == 0 for A in restricted):
            root = z
            break
    if root is None:
        raise AssertionError('Chevalley--Warning seed root was not found')
    w = h.matvec(S0, root, q)
    assert any(w) and all(h.quad(A, w, q) == 0 for A in forms)
    return S0, w, tested


def seed_disjoint_from_oil(F: GFq3, S_secret: KMatrix, S0: Matrix, V: int, q: int) -> bool:
    """Diagnostic only: test whether the public seed space meets hidden oil."""
    cols = [list(col) for col in zip(*S0)]
    projected: List[Vector] = []
    for col in cols:
        central = kmat_vec(F, S_secret, base_to_k(col))
        projected.append(k_to_base(central[:V]))
    vinegar_projection = h.columns_to_matrix(projected, 3 * V)
    return h.rank(vinegar_projection, q) == len(cols)


def construct_common_isotropic_anchored(forms: Sequence[Matrix], anchor: Vector, r: int, q: int, rng: random.Random):
    s = len(forms)
    n = len(anchor)
    W_cols: List[Vector] = [anchor[:]]
    stats = h.ConstructionStats()
    assert all(h.quad(A, anchor, q) == 0 for A in forms)
    for j in range(1, r):
        constraints: Matrix = []
        for A in forms:
            for w in W_cols:
                constraints.append(h.vecmat(w, A, q))
        L_basis = h.nullspace(constraints, n, q)
        extended = h.independent_extend(W_cols, L_basis, q)
        comp = extended[len(W_cols):]
        qdim = len(comp)
        stats.min_quotient_dim = min(stats.min_quotient_dim, qdim)
        stats.max_quotient_dim = max(stats.max_quotient_dim, qdim)
        d = 2 * s + 1
        if qdim < d:
            raise RuntimeError(f'quotient dimension {qdim} < {d}')
        # Random d-subspace of the quotient; exhaustive projective search is deterministic once chosen.
        Smat = h.random_subspace_from_complement(comp, d, q, rng)
        restricted = [h.restrict_form(A, Smat, q) for A in forms]
        root = None
        for z in h.projective_points(d, q):
            stats.projective_points_tested += 1
            if all(h.quad(A, z, q) == 0 for A in restricted):
                root = z
                break
        if root is None:
            raise AssertionError('Chevalley-Warning root was not found')
        v = h.matvec(Smat, root, q)
        W_cols.append(v)
        W = h.columns_to_matrix(W_cols, n)
        assert h.rank(W, q) == len(W_cols)
        assert all(all(x == 0 for row in h.restrict_form(A, W, q) for x in row) for A in forms)
        stats.steps += 1
    if stats.min_quotient_dim == 10**9:
        stats.min_quotient_dim = n
    return h.columns_to_matrix(W_cols, n), stats


def brute_force_projective_residual(
    zero_forms: Sequence[Matrix], scale_form: Matrix, target_scalar: int,
    U: Matrix, q: int
) -> Tuple[Vector | None, int]:
    """Solve the homogeneous residual modulo projective scaling.

    The zero_forms cut a projective zero-dimensional system in P(U).  For
    each projective root [z], the final equation asks for a nonzero scalar
    lambda with lambda^2*scale_form(z)=target_scalar.  Returning lambda*z
    therefore recovers an affine preimage without solving an extra quadratic
    equation in the projective core.
    """
    r = len(U[0]) if U and U[0] else 0
    restricted_zero = [h.restrict_form(A, U, q) for A in zero_forms]
    restricted_scale = h.restrict_form(scale_form, U, q)
    target_scalar %= q
    if target_scalar == 0:
        raise ValueError('projective residual requires a nonzero target scalar')
    tested = 0
    for z in h.projective_points(r, q):
        tested += 1
        if not all(h.quad(A, z, q) == 0 for A in restricted_zero):
            continue
        a = h.quad(restricted_scale, z, q)
        if a == 0:
            continue
        for lam in range(1, q):
            if (lam * lam * a - target_scalar) % q == 0:
                out = [(lam * zi) % q for zi in z]
                assert all(h.quad(A, out, q) == 0 for A in restricted_zero)
                assert h.quad(restricted_scale, out, q) == target_scalar
                return out, tested
    return None, tested


@dataclass
class Trial:
    success: bool
    attempts: int
    seed_subspaces: int
    final_seed_disjoint: bool
    off_oil_anchor: bool
    verified: bool
    seed_points: int
    cw_points: int
    residual_points: int


def run_trial(F: GFq3, V: int, O: int, m: int, max_attempts: int, rng: random.Random) -> Trial:
    formsK, S, forms = generate_qruov(F, V, O, m, rng)
    q = F.q
    n = 3 * (V + O)
    r = m - 3
    while True:
        y = [rng.randrange(q) for _ in range(m)]
        if any(y):
            break
    T = h.target_annihilating_transform(y, q)
    mixed = h.output_mix(forms, T, q)
    ty = h.matvec(T, y, q)
    assert all(a == 0 for a in ty[:-1]) and ty[-1] != 0

    total_seed = total_cw = total_resid = 0
    last_off = False
    last_disjoint = False
    for attempt in range(1, max_attempts + 1):
        S0, w, seed_tested = seed_anchor(mixed[:3], q, rng)
        total_seed += seed_tested
        disjoint = seed_disjoint_from_oil(F, S, S0, V, q)
        central_w = kmat_vec(F, S, base_to_k(w))
        off = any(a != F.z() for a in central_w[:V])
        if disjoint and not off:
            raise AssertionError('a disjoint seed space produced an oil anchor')
        last_disjoint = disjoint
        last_off = off
        U, stats = construct_common_isotropic_anchored(mixed[:3], w, r, q, rng)
        total_cw += stats.projective_points_tested
        z, tested = brute_force_projective_residual(
            mixed[3:-1], mixed[-1], ty[-1], U, q
        )
        total_resid += tested
        if z is None:
            continue
        x = h.matvec(U, z, q)
        ok = h.eval_map(forms, x, q) == y
        if not ok:
            raise AssertionError('candidate failed original public map')
        # Independent evaluator through K representation.
        assert public_eval(F, formsK, S, x) == y
        return Trial(True, attempt, attempt, disjoint, off, True, total_seed, total_cw, total_resid)
    return Trial(False, max_attempts, max_attempts, last_disjoint, last_off, False, total_seed, total_cw, total_resid)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, default=3)
    ap.add_argument('--V', type=int, default=3)
    ap.add_argument('--O', type=int, default=2)
    ap.add_argument('--trials', type=int, default=100)
    ap.add_argument('--max-attempts', type=int, default=10)
    ap.add_argument('--seed', type=int, default=20260729)
    ap.add_argument('--json', type=str, default=None)
    args = ap.parse_args()
    m = 3 * args.O
    n = 3 * (args.V + args.O)
    r = m - 3
    k = r - 1
    if n < 4 * r + 3:
        raise SystemExit(f'profile fails n >= 4r+3: {n} < {4*r+3}')
    F = GFq3(args.q)
    rng = random.Random(args.seed)
    results = [run_trial(F, args.V, args.O, m, args.max_attempts, rng) for _ in range(args.trials)]
    succ = sum(t.success for t in results)
    print('Exact degree-3 QR-UOV target-orthogonal toy')
    print(f'q={args.q}, V={args.V}, O={args.O}, n={n}, m={m}, isotropic_r={r}, projective_core={k}, trials={args.trials}')
    print(f'success={succ}/{args.trials}; first-attempt={sum(t.success and t.attempts==1 for t in results)}/{args.trials}')
    print(f'all returned candidates verified={all((not t.success) or t.verified for t in results)}')
    print(f'final seed spaces disjoint from oil={sum(t.final_seed_disjoint for t in results)}/{args.trials}')
    print(f'off-oil final anchors={sum(t.off_oil_anchor for t in results)}/{args.trials}')
    print(f'total seed subspaces={sum(t.seed_subspaces for t in results)}; total seed points={sum(t.seed_points for t in results)}; total extension points={sum(t.cw_points for t in results)}; total residual points={sum(t.residual_points for t in results)}')
    payload = {'parameters': vars(args) | {'m': m, 'n': n, 'isotropic_r': r, 'projective_core': k}, 'results': [asdict(t) for t in results]}
    if args.json:
        Path(args.json).write_text(json.dumps(payload, indent=2), encoding='utf-8')
        print(f'wrote {args.json}')


if __name__ == '__main__':
    main()
