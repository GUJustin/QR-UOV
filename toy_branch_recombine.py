import sympy as sp
from itertools import product
q=101
x,y,t,u=sp.symbols('x y t u')
# fixed random-ish dense quadrics over F_q
f1=3*x*x+5*x*y+7*y*y+11*x+13*y+17
f2=19*x*x+23*x*y+29*y*y+31*x+37*y+41
g1=x*x-x
g2=y*y-y
H1=sp.expand((1-t)*g1+t*f1)
H2=sp.expand((1-t)*g2+t*f2)
lam1,lam2=2,3
# eliminate y using linear form u=2x+3y
ysub=sp.expand((u-lam1*x)*pow(lam2,-1,q))
# reduce rational coefficient mod q by Poly modulus after substitution
E1=sp.Poly(sp.expand(H1.subs(y,ysub)),x,t,u,modulus=q).as_expr()
E2=sp.Poly(sp.expand(H2.subs(y,ysub)),x,t,u,modulus=q).as_expr()
R=sp.resultant(E1,E2,x)
Rpoly=sp.Poly(R,u,t,modulus=q)
# normalize monic in u over F_q(t): leading coefficient may t-dependent; at t=0 should constant unit.
# We compare relation up to unit by scale at t=0 and low t series; better extract coefficients and invert leading coeff series.
print('resultant degrees',Rpoly.degree(u),Rpoly.degree(t))

P=12
# series arrays length P mod q
def add(a,b): return [(a[i]+b[i])%q for i in range(P)]
def sub(a,b): return [(a[i]-b[i])%q for i in range(P)]
def mul(a,b):
    c=[0]*P
    for i,ai in enumerate(a):
      if ai:
       for j,bj in enumerate(b[:P-i]): c[i+j]=(c[i+j]+ai*bj)%q
    return c
def scale(a,s): return [(s*z)%q for z in a]
def inv_series(a):
    assert a[0]%q
    b=[0]*P;b[0]=pow(a[0],-1,q)
    for n in range(1,P):
      b[n]=(-b[0]*sum(a[i]*b[n-i] for i in range(1,n+1)))%q
    return b

# evaluate target quadrics on series X,Y
def evalf(X,Y,which):
    XX=mul(X,X);XY=mul(X,Y);YY=mul(Y,Y)
    if which==1:
      coeff=[3,5,7,11,13,17]
    else: coeff=[19,23,29,31,37,41]
    out=[0]*P
    for arr,c in [(XX,coeff[0]),(XY,coeff[1]),(YY,coeff[2]),(X,coeff[3]),(Y,coeff[4])]: out=add(out,scale(arr,c))
    out[0]=(out[0]+coeff[5])%q
    return out

def branch(e1,e2):
    X=[0]*P;Y=[0]*P;X[0]=e1;Y[0]=e2
    s1=(2*e1-1)%q;s2=(2*e2-1)%q
    # coefficient recurrence using H=g+t(f-g)
    for n in range(1,P):
      # recompute known series with current coeff n zero
      GX=sub(mul(X,X),X); GY=sub(mul(Y,Y),Y)
      F1=evalf(X,Y,1);F2=evalf(X,Y,2)
      known1=(GX[n]+(F1[n-1]-GX[n-1]))%q
      known2=(GY[n]+(F2[n-1]-GY[n-1]))%q
      X[n]=(-known1*pow(s1,-1,q))%q
      Y[n]=(-known2*pow(s2,-1,q))%q
    return X,Y
branches=[]
for e in product([0,1],repeat=2): branches.append(branch(*e))
# Q(U)=prod(U-u_j(t)), coefficients low-to-high in U, each series
Q=[[0]*P for _ in range(5)];Q[0][0]=1
curdeg=0
for Xs,Ys in branches:
    us=add(scale(Xs,lam1),scale(Ys,lam2))
    New=[[0]*P for _ in range(5)]
    for j in range(curdeg+1):
       New[j]=sub(New[j],mul(Q[j],us))
       New[j+1]=add(New[j+1],Q[j])
    Q=New;curdeg+=1
# resultant coefficient series in t, normalize by leading u coefficient series
coeff_expr=[]
for j in range(5):
    expr=0
    for (ju,jt),cc in Rpoly.terms():
        if ju==j:
            expr += int(cc)*t**jt
    coeff_expr.append(expr)
def expr_series(expr):
    pp=sp.Poly(expr,t,modulus=q)
    arr=[0]*P
    for (k,),c in pp.terms():
      if k<P: arr[k]=int(c)%q
    return arr
RR=[expr_series(e) for e in coeff_expr]
lead=RR[4]; invlead=inv_series(lead);RR=[mul(a,invlead) for a in RR]
for j in range(5):
    ok=Q[j]==RR[j]
    print('coef',j,ok,'Q',Q[j][:6],'R',RR[j][:6])
assert all(Q[j]==RR[j] for j in range(5))
# derivative identity on formal roots: dQ/dlambda_i at U=u_s equals -x_i Q'
# finite difference symbolically at local product level direct formula check for branch0 series.
s=0;Xs,Ys=branches[s];us=add(scale(Xs,lam1),scale(Ys,lam2))
# Q'(us)
Qp=[scale(Q[j],j) for j in range(1,5)]
def eval_poly(coeff,z):
    out=[0]*P
    for a in reversed(coeff): out=add(mul(out,z),a)
    return out
qpval=eval_poly(Qp,us)
# derivative wrt lam1 = -sum X_j prod_{l!=j}(U-u_l), evaluate at us; only j=s survives formally because factors U-u_s for j!=s.
# direct expected
expected=scale(mul(Xs,qpval),-1)
# construct derivative polynomial coefficients using product excluding each branch, then evaluate at us
Dlam=[[0]*P for _ in range(4)]
for idx,(Xj,Yj) in enumerate(branches):
  poly=[[0]*P for _ in range(4)];poly[0][0]=1;deg=0
  for ell,(Xe,Ye) in enumerate(branches):
    if ell==idx: continue
    ue=add(scale(Xe,lam1),scale(Ye,lam2));new=[[0]*P for _ in range(4)]
    for jj in range(deg+1):
      new[jj]=sub(new[jj],mul(poly[jj],ue));new[jj+1]=add(new[jj+1],poly[jj])
    poly=new;deg+=1
  for jj in range(4):Dlam[jj]=sub(Dlam[jj],mul(Xj,poly[jj]))
actual=eval_poly(Dlam,us)
print('derivative identity',actual==expected,actual[:6],expected[:6])
