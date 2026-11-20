"""Symbolic check of the quadratic tensor action in EsGB cosmology.
Action: S = int d^4x sqrt(-g) [ R/2 + f(phi) G + L_matter ],  M_pl = 1.
Metric: ds^2 = -dt^2 + a(t)^2 [ e^{h} dx^2 + e^{-h} dy^2 + dz^2 ],  h = h(t,z)  (plus polarisation).
Since det(g_ij) = a^6 exactly, matter and the minimally coupled scalar do not contribute at order h^2.
"""
import sympy as sp
from sympy.calculus.euler import euler_equations
t,z,x,y,eps=sp.symbols('t z x y epsilon')
a=sp.Function('a')(t); f=sp.Function('f')(t); h=sp.Function('h')(t,z)
X=[t,x,y,z]
g=sp.diag(-1, a**2*sp.exp(eps*h), a**2*sp.exp(-eps*h), a**2)
gi=g.inv()
n=4
Gam=[[[sum(gi[i,l]*(sp.diff(g[l,j],X[k])+sp.diff(g[l,k],X[j])-sp.diff(g[j,k],X[l])) for l in range(n))/2 for k in range(n)] for j in range(n)] for i in range(n)]
def Riem(i,j,k,l): # R^i_{jkl}
    r=sp.diff(Gam[i][j][l],X[k])-sp.diff(Gam[i][j][k],X[l])
    r+=sum(Gam[i][k][m]*Gam[m][j][l]-Gam[i][l][m]*Gam[m][j][k] for m in range(n))
    return r
Rm={}
for i in range(n):
    for j in range(n):
        for k in range(n):
            for l in range(k+1,n):
                v=sp.simplify(Riem(i,j,k,l)); Rm[(i,j,k,l)]=v; Rm[(i,j,l,k)]=-v
def R4(i,j,k,l): return Rm.get((i,j,k,l),0)
Ric=sp.Matrix(n,n,lambda j,l: sum(R4(i,j,i,l) for i in range(n)))
Rs=sum(gi[j,j]*Ric[j,j] for j in range(n))
# diagonal metric: indices are raised and lowered with the diagonal factors
Ric2=sum(gi[i,i]*gi[j,j]*Ric[i,j]**2 for i in range(n) for j in range(n))
Riem2=sum(g[i,i]*gi[j,j]*gi[k,k]*gi[l,l]*R4(i,j,k,l)**2 for i in range(n) for j in range(n) for k in range(n) for l in range(n))
GB=Rs**2-4*Ric2+Riem2
sqrtg=a**3
L=sqrtg*(Rs/2+f*GB)
L2=sp.expand(sp.diff(L,eps,2).subs(eps,0)/2)
# background check: G_FLRW = 24 H^2 (addot/a)
H=sp.diff(a,t)/a
print('background G check (must be 0):',sp.simplify(GB.subs(eps,0)-24*H**2*sp.diff(a,t,2)/a))
eq=euler_equations(L2,[h],[t,z])[0].lhs
eq=sp.expand(eq)
hd=sp.diff(h,t); hdd=sp.diff(h,t,2); hzz=sp.diff(h,z,2)
A=sp.simplify(eq.coeff(hdd)); B=sp.simplify(eq.coeff(hd)); C=sp.simplify(eq.coeff(hzz)); 
rest=sp.simplify(eq-A*hdd-B*hd-C*hzz)
GT=1+8*H*sp.diff(f,t); FT=1+8*sp.diff(f,t,2)
print('remainder (must be 0):',rest)
print('A/(a^3 G_T) =',sp.simplify(A/(a**3*GT)))
print('B - d/dt(A) =',sp.simplify(B-sp.diff(A,t)))
print('C/(a F_T)   =',sp.simplify(C/(a*FT)))
# luminal branch: fddot = H fdot  <=>  fdot = c a
c=sp.symbols('c')
print('c_T^2 - 1 on fdot = c a :',sp.simplify((FT/GT-1).subs(sp.diff(f,t,2),c*sp.diff(a,t)).subs(sp.diff(f,t),c*a)))
print('G_T on fdot = c a       :',sp.simplify(GT.subs(sp.diff(f,t),c*a)))
