"""How well do the (Xi0, n) and c_M templates capture d_GW/d_EM of the luminal Gauss-Bonnet sector?"""
import numpy as np
from scipy.optimize import least_squares
Om=0.3; OL=1-Om
E=lambda z: np.sqrt(Om*(1+z)**3+OL)
gb=lambda z,b: np.sqrt((1+b)/(1+b*E(z)/(1+z)))
xi=lambda z,X0,n: X0+(1-X0)/(1+z)**n
cm=lambda z,c: np.exp(c/(2*OL)*np.log((1+z)/(Om*(1+z)**3+OL)**(1/3)))
for zmax in (1.0,2.0,3.0):
    z=np.linspace(0.01,zmax,400)
    print(f'--- range 0 < z < {zmax}')
    for b in (0.1,0.3,0.6):
        y=gb(z,b); amp=np.max(np.abs(y-1))
        r1=least_squares(lambda p: xi(z,*p)-y,[1.01,2.0],bounds=([0.5,0.01],[1.5,20]))
        r2=least_squares(lambda p: cm(z,p[0])-y,[0.0])
        print(f' beta={b}: max anomaly {amp*100:.2f}% | (Xi0,n)=({r1.x[0]:.4f},{r1.x[1]:.2f}) max residual {np.max(np.abs(r1.fun))*100:.2f}% | c_M={r2.x[0]:+.3f} max residual {np.max(np.abs(r2.fun))*100:.2f}%')
