#!/usr/bin/env python3
"""Toy check for specialization from a generic to a colliding separator.

Four explicit polynomial branches are combined by product/interpolation formulas.
The separator x+eps*y is generic, while its specialization x has a collision at
t=0.  The symmetric products specialize coefficientwise and derivative recovery
at t=1 remains exact for every simple target root.
"""
import sympy as sp

t,U,e=sp.symbols('t U e')
branches=[(t,0),(-t,1),(t+2,3),(t+3,4)]

def build(eps):
    us=[sp.expand(x+eps*y) for x,y in branches]
    bs=[sp.expand(y+2*x) for x,y in branches]
    q=sp.expand(sp.prod(U-u for u in us))
    v=sp.expand(sum(b*sp.prod(U-us[j] for j in range(len(us)) if j!=i)
                    for i,b in enumerate(bs)))
    return us,bs,q,v

ug,bg,qg,vg=build(e)
u0,b0,q0,v0=build(0)
assert sp.expand(qg.subs(e,0)-q0)==0
assert sp.expand(vg.subs(e,0)-v0)==0
assert sp.expand(u0[0].subs(t,0)-u0[1].subs(t,0))==0  # start collision
qt=sp.Poly(q0.subs(t,1),U)
assert sp.gcd(qt,qt.diff()).degree()==0
for u,b in zip(u0,b0):
    tau=sp.expand(u.subs(t,1)); want=sp.expand(b.subs(t,1))
    got=sp.cancel(v0.subs({t:1,U:tau})/sp.diff(q0,U).subs({t:1,U:tau}))
    assert sp.expand(got-want)==0
print('generic-to-local specialization: PASS')
print('start separator values:',[x.subs(t,0) for x in u0])
print('target separator values:',[x.subs(t,1) for x in u0])
print('deg_t q_generic / q_local:',sp.Poly(qg,t).degree(),sp.Poly(q0,t).degree())
print('deg_t v_generic / v_local:',sp.Poly(vg,t).degree(),sp.Poly(v0,t).degree())
print('simple target derivative recovery: PASS')
