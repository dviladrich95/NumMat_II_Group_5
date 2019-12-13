import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt

def u_klm(kk,ll,mm,xx,yy,zz):
    a_klm=np.sqrt(8)

    lambda_klm=np.pi**2*(kk**2+ll**2+mm**2)
    z_klm=8*a_klm**3/(np.pi**3*kk*ll*mm*lambda_klm)
    return z_klm*np.sin(np.pi*kk*xx)*np.sin(np.pi*ll*yy)*np.sin(np.pi*mm*zz)

N=30



R=10




def a07ex02getsola(Nx,Ny,Nz,R):


    x = np.linspace(0, 1, Nx)
    y = np.linspace(0, 1, Ny)
    z = np.linspace(0, 1, Nz)

    k=2*np.linspace(0,R//2,R//2+1)+1
    l=2*np.linspace(0,R//2,R//2+1)+1
    m=2*np.linspace(0,R//2,R//2+1)+1

    kk,ll,mm,xx,yy,zz=np.meshgrid(k,l,m,x,y,z)

    u_klm_func=u_klm(kk,ll,mm,xx,yy,zz)

    u_h=np.sum(u_klm_func,(0,1,2))

    return u_h


#sns.heatmap(slice)
#plt.show()