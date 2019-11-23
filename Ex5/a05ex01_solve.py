import numpy as np
import matplotlib.pyplot as plt

from a05ex01_get_laplace import a05ex01_get_laplace

N =  32
L = 2.2


def is_in_domain(xx,yy,L):
    return (xx-L/2)**2 + (yy-L/2)**2 <= 1
def f(xx,yy):
    return 0*xx+1
def g(xx,yy):
    return 0*xx

    
Lh,xx,yy,bar_omega_h,omega_h,gamma_h = a05ex01_get_laplace(L,N,is_in_domain)

plt.figure()
plt.plot(xx[np.logical_not(bar_omega_h)],yy[np.logical_not(bar_omega_h)],'yv')
plt.plot(xx[bar_omega_h],yy[bar_omega_h],'ko')
plt.show()


# u   = np.zeros((N,N)) # in MATLAB setting nan allows to show only domain
# rhs = np.zeros((N,N))
# rhs[omega_h]   = +f(xx[omega_h],yy[omega_h]);
# rhs[gamma_h]   = +g(xx[gamma_h],yy[gamma_h]);
# u[bar_omega_h] = sp.linalg.spsolve(Lh, rhs[bar_omega_h])