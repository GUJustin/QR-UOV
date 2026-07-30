#!/usr/bin/env python3
"""Exact toy for the corrected infinite-branch specialization argument.

One branch u_inf=1/t escapes, while u=2 and u=3 stay finite. After multiplying
q and its scalar interpolation numerator by the common power t, the escaping
branch's numerator summand does not vanish. It is, however, divisible by the
entire finite eliminant, so evaluation at either finite simple root recovers the
correct scalar value.
"""
import sympy as sp

t, U = sp.symbols("t U")
branches = [1 / t, sp.Integer(2), sp.Integer(3)]
values = list(branches)
q = sp.prod(U - u for u in branches)
v = sum(b * sp.prod(U - branches[j] for j in range(len(branches)) if j != i)
        for i, b in enumerate(values))
qbar = sp.expand(sp.limit(t * q, t, 0))
vbar = sp.expand(sp.limit(t * v, t, 0))
infinite_term = sp.expand((U - 2) * (U - 3))
assert sp.rem(infinite_term, qbar, domain=sp.QQ) == 0
for tau in (2, 3):
    recovered = sp.simplify(vbar.subs(U, tau) / sp.diff(qbar, U).subs(U, tau))
    assert recovered == tau, (tau, recovered)
print("Infinite-branch specialization toy: PASS")
print("  normalized eliminant:", qbar)
print("  normalized numerator:", vbar)
print("  escaping-branch summand:", infinite_term, "(divisible by eliminant)")
print("  recovered finite values: 2, 3")
