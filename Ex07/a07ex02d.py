import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt


def w_klm(kk,ll,mm):
    a_klm=np.sqrt(8)
    lambda_klm=np.pi**2*(kk**2+ll**2+mm**2)
    w_klm=8*a_klm/(np.pi**3*kk*ll*mm*lambda_klm)
    return w_klm

def wt_klm(T,kk,ll,mm):
    lambda_klm=np.pi**2*(kk**2+ll**2+mm**2)
    t_fact=1 - np.exp(-lambda_klm * T)
    w=w_klm(kk, ll, mm)
    return w*t_fact

def u_klm(kk,ll,mm,xx,yy,zz):
    a_klm=np.sqrt(8)
    return a_klm*np.sin(np.pi*kk*xx)*np.sin(np.pi*ll*yy)*np.sin(np.pi*mm*zz)

def a07ex02getsold(T,Nx,Ny,Nz,R):

    x = np.linspace(0, 1, Nx+2)
    y = np.linspace(0, 1, Ny+2)
    z = np.linspace(0, 1, Nz+2)

    k=2*np.linspace(0,R//2,R//2+1)+1
    l=2*np.linspace(0,R//2,R//2+1)+1
    m=2*np.linspace(0,R//2,R//2+1)+1

    kk,ll,mm,xx,yy,zz=np.meshgrid(k,l,m,x,y,z)

    wt=wt_klm(T, kk, ll, mm)
    u_h=np.sum(wt*u_klm(kk,ll,mm,xx,yy,zz),(0,1,2))

    return u_h