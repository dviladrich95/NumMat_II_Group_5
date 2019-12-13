import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt

def u_klm(kk,ll,xx,yy):#,zz):
    a_klm=np.sqrt(4)

    lambda_klm=np.pi**2*(kk**2+ll**2)#+mm**2)
    z_klm=8*a_klm/(np.pi**2*kk*ll*lambda_klm)
    return a_klm**2*np.sin(np.pi*kk*xx)**2*np.sin(np.pi*ll*yy)**2#*np.sin(np.pi*mm*zz)**2

N=200

Nx=N
Ny=N
Nz=N

R=5

x=np.linspace(0,1,Nx)
y=np.linspace(0,1,Ny)
#z=np.linspace(0,1,Nz)

k=2*np.linspace(0,R//2,R//2)+1
l=2*np.linspace(0,R//2,R//2)+1
#m=2*np.linspace(0,R//2,R//2)+1

kk,ll,xx,yy=np.meshgrid(k,l,x,y)

print(kk.shape,xx.shape)
u_klm_func=u_klm(kk,ll,xx,yy)#,zz)


for i,dim_1_slice in enumerate(u_klm_func):
    print(i)
    for k,dim_3_slice in enumerate(dim_1_slice):
        print(i,k)
        print(dim_3_slice.shape)
        print(np.sum(dim_3_slice)/200**2)
        sns.heatmap(dim_3_slice)
        plt.show()

u_h=np.sum(u_klm_func,(0,1))

sns.heatmap(u_h)
plt.show()