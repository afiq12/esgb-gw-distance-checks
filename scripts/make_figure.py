"""Figure: d_GW/d_EM for the luminal Gauss-Bonnet sector versus the non-minimally coupled model."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
Om=0.3; z=np.linspace(0,3,400); aH=np.sqrt(Om*(1+z)**3+1-Om)/(1+z)
gb=lambda b: np.sqrt((1+b)/(1+b*aH)); xi=lambda X0,n: X0+(1-X0)/(1+z)**n
fig,ax=plt.subplots(figsize=(6.4,3.3),dpi=220)
ax.fill_between(z,xi(0.8,3.8),xi(1.3,3.8),color='#bbbbbb',alpha=.35,lw=0,label=r'GWTC-5.0 interval quoted by Lagos & Wolf ($\Xi_0=1.0^{+0.3}_{-0.2}$, narrow $H_0$ prior; $n=3.8$)')
ax.plot(z,xi(0.88,3.2),color='#c0392b',lw=1.8,label=r'non-minimal coupling (Lagos & Wolf 2026): $\Xi_0=0.88$, $n=3.2$')
for b,ls in [(0.1,':'),(0.3,'--'),(0.6,'-')]: ax.plot(z,gb(b),color='#1f4e79',lw=1.8,ls=ls,label=rf'luminal Gauss-Bonnet, $\beta={b}$')
ax.axhline(1,color='k',lw=.5); ax.set_xlabel('redshift $z$'); ax.set_ylabel(r'$d_{\rm GW}/d_{\rm EM}$'); ax.set_xlim(0,3); ax.set_ylim(0.78,1.32)
ax.legend(fontsize=6.5,loc='upper left',frameon=False); fig.tight_layout(); fig.savefig('ratio_gb_vs_nmc.png')
