import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt

def u_klm(kk,ll,mm,xx,yy,zz):
    a_klm=np.sqrt(8)

    lambda_klm=np.pi**2*(kk**2+ll**2+mm**2)
    z_klm=8*a_klm/(np.pi**3*kk*ll*mm*lambda_klm)
    return z_klm*np.sin(np.pi*kk*xx)*np.sin(np.pi*ll*yy)*np.sin(np.pi*mm*zz)

N=70

Nx=N
Ny=N
Nz=N

R=10

x=np.linspace(0,1,Nx)
y=np.linspace(0,1,Ny)
z=np.linspace(0,1,Nz)

k=2*np.linspace(0,R//2,R//2+1)+1
l=2*np.linspace(0,R//2,R//2+1)+1
m=2*np.linspace(0,R//2,R//2+1)+1

kk,ll,mm,xx,yy,zz=np.meshgrid(k,l,m,x,y,z)

print(kk.shape,xx.shape)
u_klm_func=u_klm(kk,ll,mm,xx,yy,zz)

'''
for i,dim_1_slice in enumerate(u_klm_func):
    print(i)
    for j,dim_2_slice in enumerate(dim_1_slice):
        print(i,j)
        for k,dim_3_slice in enumerate(dim_2_slice):
            print(i,j,k)
            print(np.sum(dim_3_slice)/N**3)
            for slice in dim_3_slice:
                sns.heatmap(slice)
                plt.show()
'''

u_h=np.sum(u_klm_func,(0,1,2))

for slice in u_h:
    sns.heatmap(slice)
    plt.show()