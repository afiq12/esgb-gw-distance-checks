"""Numerical check: the GW amplitude on the luminal EsGB branch satisfies
d_GW/d_EM = sqrt(G_T(0)/G_T(z)),  G_T = 1 + beta * a H / H0,  beta = 8 c H0.
LambdaCDM background (test approximation), units H0 = 1."""
import numpy as np
from scipy.integrate import solve_ivp
Om=0.3
Hf=lambda a: np.sqrt(Om/a**3+1-Om)
def run(beta,ze,k=400.0):
    GT=lambda a: 1+beta*a*Hf(a)
    ae=1/(1+ze)
    def rhs(eta,y):
        a,h,p=y; return [a*a*Hf(a), p/(a*a*GT(a)), -a*a*GT(a)*k*k*h]
    ev=lambda eta,y: y[0]-1.0
    ev.terminal=True
    s=solve_ivp(rhs,[0,10],[ae,1.0,0.0],events=ev,rtol=1e-10,atol=1e-12,max_step=0.002)
    a,h,p=s.y[:,-1]; hp=p/(a*a*GT(a))
    return np.sqrt(h*h+(hp/k)**2)   # envelope at a=1, initial envelope = 1
print(' beta   z_e   numerical  analytic   relative difference')
for beta in (0.1,0.3,0.6):
    for ze in (0.3,0.67,1.5,3.0):
        num=run(0.0,ze)/run(beta,ze)          # d_GW/d_EM = A_GR/A_MG
        ana=np.sqrt((1+beta)/(1+beta*Hf(1/(1+ze))/(1+ze)))
        print(f' {beta:4.1f}  {ze:4.2f}  {num:.6f}  {ana:.6f}  {abs(num/ana-1):.1e}')
