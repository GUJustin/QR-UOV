#!/usr/bin/env python3
"""Toy check of split-before-lift with the Vandermonde product start.

Over F_101 in two variables, independently lift all four simple start roots,
recombine q and interpolation numerators, and compare q with direct elimination.
Also checks scalar linear/quadratic value numerators and coordinate recovery.
"""
import sympy as sp
from itertools import product
q=101
x,y,t,u=sp.symbols('x y t u')
P=14
# Vandermonde start L_a=x+a y+a^2, g1=L1L2, g2=L3L4.
def L(a): return x+a*y+a*a
g1=sp.expand(L(1)*L(2)); g2=sp.expand(L(3)*L(4))
# Dense target quadrics.
f1=3*x*x+5*x*y+7*y*y+11*x+13*y+17
f2=19*x*x+23*x*y+29*y*y+31*x+37*y+41
H1=sp.expand((1-t)*g1+t*f1); H2=sp.expand((1-t)*g2+t*f2)
lam=(2,3)
# Direct eliminant under u=2x+3y.
ysub=sp.expand((u-lam[0]*x)*pow(lam[1],-1,q))
E1=sp.Poly(H1.subs(y,ysub),x,t,u,modulus=q).as_expr()
E2=sp.Poly(H2.subs(y,ysub),x,t,u,modulus=q).as_expr()
Rpoly=sp.Poly(sp.resultant(E1,E2,x),u,t,modulus=q)
print('direct degrees',Rpoly.degree(u),Rpoly.degree(t))

def zero(): return [0]*P
def add(a,b): return [(a[i]+b[i])%q for i in range(P)]
def sub(a,b): return [(a[i]-b[i])%q for i in range(P)]
def scale(a,c): return [(c*z)%q for z in a]
def mul(a,b):
    c=[0]*P
    for i,ai in enumerate(a):
        if ai:
            for j in range(P-i): c[i+j]=(c[i+j]+ai*b[j])%q
    return c
def invs(a):
    assert a[0]%q
    b=[0]*P;b[0]=pow(a[0],-1,q)
    for n in range(1,P): b[n]=(-b[0]*sum(a[i]*b[n-i] for i in range(1,n+1)))%q
    return b

def expr_series(expr,X,Y):
    # Quadratic expression evaluator on series.
    poly=sp.Poly(expr,x,y,modulus=q)
    out=zero()
    for (i,j),cc in poly.terms():
        term=[1]+[0]*(P-1)
        for _ in range(i): term=mul(term,X)
        for _ in range(j): term=mul(term,Y)
        out=add(out,scale(term,int(cc)))
    return out

def inv2(A):
    a,b=A[0];c,d=A[1];det=(a*d-b*c)%q;di=pow(det,-1,q)
    return [[d*di%q,-b*di%q],[-c*di%q,a*di%q]]

def branch(root):
    X=[root[0]]+[0]*(P-1);Y=[root[1]]+[0]*(P-1)
    J=sp.Matrix([[sp.diff(g1,x).subs({x:root[0],y:root[1]}),sp.diff(g1,y).subs({x:root[0],y:root[1]})],
                 [sp.diff(g2,x).subs({x:root[0],y:root[1]}),sp.diff(g2,y).subs({x:root[0],y:root[1]})]])
    Ji=inv2([[int(J[i,j])%q for j in range(2)] for i in range(2)])
    for n in range(1,P):
        G=[expr_series(g1,X,Y),expr_series(g2,X,Y)]
        F=[expr_series(f1,X,Y),expr_series(f2,X,Y)]
        known=[(G[i][n]+F[i][n-1]-G[i][n-1])%q for i in range(2)]
        X[n]=(-Ji[0][0]*known[0]-Ji[0][1]*known[1])%q
        Y[n]=(-Ji[1][0]*known[0]-Ji[1][1]*known[1])%q
    return X,Y

# Choose one root of each factor pair; p(T)=T^2+yT+x.
roots=[]
for a,b in product([1,2],[3,4]): roots.append((a*b%q,(-(a+b))%q))
branches=[branch(r) for r in roots]
# Product tree, carrying coordinate and scalar value numerators.
hquad=7*x*x+9*x*y+12*y*y+5*x+4*y+3
families=5 # q plus x,y,linear,quadratic numerators
Q=[[0]*P for _ in range(5)];Q[0][0]=1
WX=[[0]*P for _ in range(5)];WY=[[0]*P for _ in range(5)]
WL=[[0]*P for _ in range(5)];WH=[[0]*P for _ in range(5)]
deg=0
for X,Y in branches:
    us=add(scale(X,lam[0]),scale(Y,lam[1]))
    vals=[X,Y,add(scale(X,6),scale(Y,8)),expr_series(hquad,X,Y)]
    newQ=[[0]*P for _ in range(5)]
    newWs=[[[0]*P for _ in range(5)] for __ in range(4)]
    oldWs=[WX,WY,WL,WH]
    for j in range(deg+1):
        newQ[j]=sub(newQ[j],mul(Q[j],us));newQ[j+1]=add(newQ[j+1],Q[j])
        for z in range(4):
            newWs[z][j]=sub(newWs[z][j],mul(oldWs[z][j],us))
            newWs[z][j]=add(newWs[z][j],mul(vals[z],Q[j]))
            newWs[z][j+1]=add(newWs[z][j+1],oldWs[z][j])
    Q=newQ;WX,WY,WL,WH=newWs;deg+=1
# Compare q to direct resultant after monic normalization.
def expr_to_series(expr):
    pp=sp.Poly(expr,t,modulus=q);a=[0]*P
    for (k,),c in pp.terms():
        if k<P:a[k]=int(c)%q
    return a
RR=[]
for j in range(5):
    e=0
    for (ju,jt),cc in Rpoly.terms():
        if ju==j:e+=int(cc)*t**jt
    RR.append(expr_to_series(e))
leadinv=invs(RR[4]);RR=[mul(a,leadinv) for a in RR]
assert Q==RR
print('eliminant match: True')
# At each formal branch u_s, W_value/Q' equals the desired value.
def peval(coeff,z):
    out=zero()
    for a in reversed(coeff):out=add(mul(out,z),a)
    return out
Qp=[scale(Q[j],j) for j in range(1,5)]
for idx,(X,Y) in enumerate(branches):
    us=add(scale(X,lam[0]),scale(Y,lam[1]));den=peval(Qp,us);dinv=invs(den)
    expected=[X,Y,add(scale(X,6),scale(Y,8)),expr_series(hquad,X,Y)]
    got=[mul(peval(W,us),dinv) for W in [WX,WY,WL,WH]]
    assert got==expected
print('coordinate/linear/quadratic numerator identities: True')
print('all checks passed')
