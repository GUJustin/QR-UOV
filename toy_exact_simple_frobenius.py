#!/usr/bin/env python3
"""Exhaustive F_9/F_3 toy for exact simple-root Frobenius descent.

F_9 = F_3[u]/(u^2+1), so u^2=2 and z^3 is conjugation.
For a Frobenius-stable target pair P,P^3, the local test
    (lambda^(3^{-1}) . P)^3 = lambda . P
passes for non-base P exactly when lambda collides on P and P^3.  Hence every
passing non-base point is attached to a repeated, not simple, separator root.
"""
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class F9:
    a: int
    b: int
    def __post_init__(self):
        object.__setattr__(self, 'a', self.a % 3)
        object.__setattr__(self, 'b', self.b % 3)
    def __add__(self, other):
        other = coerce(other); return F9(self.a+other.a, self.b+other.b)
    __radd__ = __add__
    def __neg__(self): return F9(-self.a,-self.b)
    def __sub__(self, other): return self + (-coerce(other))
    def __rsub__(self, other): return coerce(other)-self
    def __mul__(self, other):
        other=coerce(other)
        # u^2=-1=2 mod 3
        return F9(self.a*other.a+2*self.b*other.b,
                  self.a*other.b+self.b*other.a)
    __rmul__=__mul__
    def __pow__(self,n):
        out=F9(1,0); x=self
        while n:
            if n&1: out=out*x
            x=x*x; n//=2
        return out
    def __repr__(self): return f'({self.a}+{self.b}u)'

def coerce(x): return x if isinstance(x,F9) else F9(x,0)

def dot(lam,p): return sum((x*y for x,y in zip(lam,p)),F9(0,0))
def frob_vec(p): return tuple(x**3 for x in p)
def test(lam,p):
    # d=2, so q^{-1}=q on F_9.
    lam_inv_frob=frob_vec(lam)
    return dot(lam_inv_frob,p)**3 == dot(lam,p)

els=[F9(a,b) for a in range(3) for b in range(3)]
u=F9(0,1)
P=(u,F9(1,0)); Pq=frob_vec(P)
assert Pq != P and frob_vec(Pq)==P
base=(F9(1,0),F9(2,0))

passing=0; colliding=0; simple_rejects=0
for lam in ((x,y) for x in els for y in els):
    collide=dot(lam,P)==dot(lam,Pq)
    passes=test(lam,P)
    assert passes==collide
    passing += passes
    colliding += collide
    if not collide:
        assert not passes
        simple_rejects += 1
    assert test(lam,base)

# One K'-linear hyperplane among 9^2 choices has 9 elements.
assert passing==colliding==9
assert simple_rejects==72
print('exact simple-root Frobenius toy: PASS')
print('non-base target orbit:',P,Pq)
print('passing/colliding lambda:',passing,'of',81)
print('simple separator choices rejecting non-base point:',simple_rejects)
print('all lambda accept a base-field point: PASS')
