"""Preliminary estimate: how many Einstein Telescope bright sirens separate a
modified d_GW/d_EM from general relativity.

One-parameter Fisher estimate. The signal is A * ln(d_GW/d_EM)(z) with A = 1 the
model and A = 0 general relativity. Two cases: background cosmology known from
electromagnetic data, and H0 free with no prior (a constant offset in ln d is
marginalised). The truth lies between the two.

Simple error model (assumptions of this script): relative distance error
sigma^2 = sigma_inst^2 + sigma_lens^2 with sigma_inst = 0.1449 z - 0.0118 z^2 + 0.0012 z^3
and sigma_lens = 0.05 z; sirens with counterparts distributed as
P(z) ~ 4 pi d_C^2 R(z) / (H (1+z)), R = 1+2z (z<=1), 3(5-z)/4 (z>1), for 0.07 < z < 2.
These are the fitting forms commonly used in the standard-siren forecasting
literature for the Einstein Telescope; they are not derived here. Peculiar
velocities, selection effects and population uncertainties are ignored.

Models: non-minimal coupling template Xi_0 = 0.88, n = 3.2 (Lagos and Wolf,
arXiv:2609.04112) and the luminal Gauss-Bonnet ratio sqrt((1+beta)/(1+beta aH/H0)).
LambdaCDM background (Omega_m = 0.3) as a test approximation."""
import numpy as np
from scipy.integrate import quad
Om=0.3
E=lambda z: np.sqrt(Om*(1+z)**3+1-Om)
aH=lambda z: E(z)/(1+z)
gb=lambda z,b: np.sqrt((1+b)/(1+b*aH(z)))
sig=lambda z: np.sqrt((0.1449*z-0.0118*z**2+0.0012*z**3)**2+(0.05*z)**2)
dc=lambda z: quad(lambda x:1/E(x),0,z)[0]
R=lambda z: 1+2*z if z<=1 else 0.75*(5-z)
P=lambda z: 4*np.pi*dc(z)**2*R(z)/(E(z)*(1+z))
zmin,zmax=0.07,2.0
norm=quad(P,zmin,zmax)[0]
avg=lambda f: quad(lambda z:f(z)*P(z)/norm,zmin,zmax,limit=200)[0]
def info(s):
    FAA=avg(lambda z:(s(z)/sig(z))**2); FAh=avg(lambda z:s(z)/sig(z)**2); Fhh=avg(lambda z:1/sig(z)**2)
    return FAA, FAA-FAh**2/Fhh
print('model                    N(2s) N(3s) N(5s)  [cosmology known]   N(2s) N(3s) N(5s)  [H0 free]')
models=[('non-minimal Xi0=0.88 n=3.2',lambda z: np.log(0.88+0.12/(1+z)**3.2))]
models+=[(f'Gauss-Bonnet beta={b:g}',(lambda b:(lambda z:np.log(gb(z,b))))(b)) for b in (0.3,0.1,0.03)]
for lab,s in models:
    a,m=info(s)
    print(f'{lab:26s} {4/a:8.0f} {9/a:8.0f} {25/a:8.0f}      {4/m:8.0f} {9/m:8.0f} {25/m:8.0f}')
kpk=0.0664
a,m=info(lambda z:0.5*(1-aH(z))/kpk)
print(f'Gauss-Bonnet, small amplitude: N(2 sigma) = {4/a:.3f}/amp^2 (cosmology known), {4/m:.3f}/amp^2 (H0 free)')
for lab,amp in [('BBN ceiling, T=0.1 MeV',1.09e-9),('recombination, assumed 5% tolerance',1.59e-4)]:
    print(f'  {lab}: amp = {amp:.2e}  ->  N(2 sigma) between {4/a/amp**2:.1e} and {4/m/amp**2:.1e}')
