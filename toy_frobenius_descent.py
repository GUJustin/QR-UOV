#!/usr/bin/env python3
"""Exhaustive toy check of the two-scalar Frobenius descent test.

Works in F_25 = F_5[a]/(a^2+2).  For a non-F_5 point P, a uniform
rho in F_25^2 should satisfy rho.(P^5-P)=0 with probability exactly 1/25.
The test is evaluated using only the two scalar values rho.P and rho^(5).P.
"""
from __future__ import annotations
from dataclasses import dataclass
from itertools import product

Q = 5

@dataclass(frozen=True)
class F25:
    a: int
    b: int
    def __post_init__(self):
        object.__setattr__(self, 'a', self.a % Q)
        object.__setattr__(self, 'b', self.b % Q)
    def __add__(self, other): return F25(self.a+other.a, self.b+other.b)
    def __sub__(self, other): return F25(self.a-other.a, self.b-other.b)
    def __neg__(self): return F25(-self.a,-self.b)
    def __mul__(self, other):
        # alpha^2 = -2 = 3 mod 5
        return F25(self.a*other.a + 3*self.b*other.b,
                   self.a*other.b + self.b*other.a)
    def __pow__(self, n: int):
        out=F25(1,0); base=self
        while n:
            if n&1: out=out*base
            base=base*base; n//=2
        return out

ZERO=F25(0,0)
ELTS=[F25(a,b) for a in range(Q) for b in range(Q)]

def dot(x,y):
    out=ZERO
    for u,v in zip(x,y): out=out+u*v
    return out

def frob(x): return x**Q

P_nonbase=(F25(0,1),F25(1,1))
P_base=(F25(2,0),F25(4,0))

for name,P in [('nonbase',P_nonbase),('base',P_base)]:
    passed=0
    total=0
    for rho in product(ELTS, repeat=2):
        # d=2, so q^{-1}=q on F_25.
        left=frob(dot(tuple(frob(r) for r in rho),P))
        right=dot(rho,P)
        passed += (left==right)
        total += 1
    print(name, passed, total, passed/total)

assert sum(
    frob(dot(tuple(frob(r) for r in rho),P_nonbase))==dot(rho,P_nonbase)
    for rho in product(ELTS, repeat=2)
)==Q**2
assert all(
    frob(dot(tuple(frob(r) for r in rho),P_base))==dot(rho,P_base)
    for rho in product(ELTS, repeat=2)
)
print('all checks passed')
