"""(1) EsGB Friedmann equation from the minisuperspace action with lapse N.
(2) Check of G_T and F_T through the Horndeski map (Kobayashi, Yamaguchi, Yokoyama 2011)."""
import sympy as sp
t=sp.symbols('t'); a=sp.Function('a')(t); N=sp.Function('N')(t); f=sp.Function('f')(t); ph=sp.Function('phi')(t)
x,y,z=sp.symbols('x y z'); X=[t,x,y,z]; n=4
g=sp.diag(-N**2,a**2,a**2,a**2); gi=g.inv()
Gam=[[[sum(gi[i,l]*(sp.diff(g[l,j],X[k])+sp.diff(g[l,k],X[j])-sp.diff(g[j,k],X[l])) for l in range(n))/2 for k in range(n)] for j in range(n)] for i in range(n)]
def Riem(i,j,k,l):
    return sp.simplify(sp.diff(Gam[i][j][l],X[k])-sp.diff(Gam[i][j][k],X[l])+sum(Gam[i][k][m]*Gam[m][j][l]-Gam[i][l][m]*Gam[m][j][k] for m in range(n)))
R4={(i,j,k,l):Riem(i,j,k,l) for i in range(n) for j in range(n) for k in range(n) for l in range(n)}
Ric=sp.Matrix(n,n,lambda j,l: sum(R4[(i,j,i,l)] for i in range(n)))
Rs=sum(gi[j,j]*Ric[j,j] for j in range(n))
Ric2=sum(gi[i,i]*gi[j,j]*Ric[i,j]**2 for i in range(n) for j in range(n))
Riem2=sum(g[i,i]*gi[j,j]*gi[k,k]*gi[l,l]*R4[(i,j,k,l)]**2 for i in range(n) for j in range(n) for k in range(n) for l in range(n))
GB=sp.simplify(Rs**2-4*Ric2+Riem2); sq=N*a**3
print('sqrt(-g) G - d/dt(8 adot^3/N^3) =',sp.simplify(sq*GB-sp.diff(8*sp.diff(a,t)**3/N**3,t)))
V=sp.Function('V')(ph); rho_m=sp.symbols('rho_m')
# minisuperspace Lagrangian (after integration by parts); matter enters as -N a^3 rho_m
Lm=-3*a*sp.diff(a,t)**2/N-8*sp.diff(f,t)*sp.diff(a,t)**3/N**3+a**3*(sp.diff(ph,t)**2/(2*N)-N*V)-N*a**3*rho_m
con=sp.simplify(sp.diff(Lm,N).subs(N,1).doit()/a**3)
H=sp.symbols('H'); con=sp.simplify(con.subs(sp.diff(a,t),H*a))
print('Hamiltonian constraint / a^3 =',sp.expand(con),' = 0')
# (2) Horndeski map
Xk,lnX,Hh,pd,pdd=sp.symbols('X lnX H phidot phiddot'); f1,f2=sp.symbols('f1 f2') # f', f''
G4=sp.Rational(1,2)+4*f2*Xk*(2-lnX); G4X=4*f2*(2-lnX)+4*f2*Xk*(-1/Xk)
G5X=-4*f1/Xk; G5p=-4*f2*lnX
GT=2*(G4-2*Xk*G4X-Xk*(Hh*pd*G5X-G5p)); FT=2*(G4-Xk*(pdd*G5X+G5p))
print('G_T =',sp.simplify(GT),'  (= 1 + 8 H fdot, fdot = f1 phidot)')
print('F_T =',sp.simplify(FT.subs(Xk,pd**2/2)),'  (= 1 + 8 fddot, fddot = f2 phidot^2 + f1 phiddot)')
