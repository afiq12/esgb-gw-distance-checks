"""Preliminary estimate: cosmology-only ceiling on the luminal Gauss-Bonnet coupling.

On the luminal branch G_T = 1 + beta * a H / H0 for any scalar potential, and the
Hamiltonian constraint reads 3 H^2 G_T = rho_tot. If the scalar energy is subdominant
at nucleosynthesis, the expansion rate is H^2 = H^2_std / G_T, which maps onto an
effective number of neutrino species through H^2/H^2_std = 1 + 7 (N_nu - 3)/43.

Input from the literature: N_nu = 2.898 +/- 0.141 (BBN+CMB; Yeh, Shelton, Olive and
Fields, JCAP 10 (2022) 046). The two-sigma lower edge is taken as N_nu - 2 sigma,
assuming a Gaussian error (an assumption of this script, not a number quoted there).

Limitations: G_T - 1 scales like the temperature during BBN, whereas the N_nu mapping
assumes a constant rescaling of H. The script therefore brackets the result between
T = 1 MeV (weak freeze-out) and T = 0.1 MeV (onset of nucleosynthesis). The
recombination line uses an ASSUMED tolerance on H^2, not a literature bound.
LambdaCDM background (Omega_m = 0.3, h = 0.674) as a test approximation."""
import numpy as np
from scipy.optimize import minimize_scalar
Om=0.3; h=0.674
E=lambda z: np.sqrt(Om*(1+z)**3+1-Om)
aH=lambda z: E(z)/(1+z)
r=minimize_scalar(lambda z:-(0.5*(1-aH(z))),bounds=(0.01,2),method='bounded')
zpk,kpk=r.x,-r.fun
print(f'small-beta anomaly: max(d_GW/d_EM - 1) = {kpk:.4f} * beta at z = {zpk:.3f}')
T0=2.7255*8.617333e-11   # CMB temperature today, MeV
Mpl=1.220890e22          # Planck mass, MeV
hbar=6.582119e-22        # MeV s
H0=h*100/3.0857e19       # 1/s
def aH_over_H0(T,g,gs):
    a=(T0/T)*(3.909/gs)**(1/3)
    H=1.6602*np.sqrt(g)*T*T/Mpl/hbar
    return a*H/H0
Nnu,sig=2.898,0.141
S2=1+7*((Nnu-2*sig)-3.0)/43
eps=1/S2-1
print(f'BBN: N_nu >= {Nnu-2*sig:.3f} (2 sigma)  ->  H^2/H^2_std >= {S2:.4f}  ->  G_T - 1 <= {eps:.4f}')
for T,g,gs in [(1.0,10.75,10.75),(0.1,3.36,3.909)]:
    x=aH_over_H0(T,g,gs); b=eps/x
    print(f'  T = {T:3.1f} MeV: aH/H0 = {x:.3e}   beta <= {b:.2e}   max anomaly <= {kpk*b:.2e}')
zr=1090; Or=4.18e-5/h**2
xr=np.sqrt(Om*(1+zr)**3+Or*(1+zr)**4)/(1+zr)
for e in (0.05,0.02):
    print(f'recombination (assumed tolerance {e}): aH/H0 = {xr:.1f}   beta <= {e/xr:.2e}   max anomaly <= {kpk*e/xr:.2e}')
