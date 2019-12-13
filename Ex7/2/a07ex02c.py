import numpy as np
from scipy.sparse.linalg import spsolve
from scipy import sparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import seaborn as sns

def a07ex02c(N,M,T):
    # parameters
    L    = 1.0
    k    = N + 2    # define k = N+2 for ease of typing
    l    = M + 2
    tau  = T/(l-1)
    h    = L/(k-1)  # lattice spacing

    # array construction
    xh, yh , zh = np.meshgrid(np.linspace(0,L,k), np.linspace(0,L,k),np.linspace(0,L,k)) # create mesh using MATLABs meshgrid

    xh = np.transpose(xh)
    yh = np.transpose(yh)
    zh = np.transpose(zh)


    ix = np.reshape(np.arange(0,k**3),(k,k,k)).T # create indices in lexic. ordering

    ixyz  = ix[1:-1,1:-1,1:-1].flatten()   # index (i,j,k)

    ixmyz = ix[0:-2,1:-1,1:-1].flatten()  # index (i-1,j,k)
    ixpyz = ix[2:,1:-1,1:-1].flatten()     # index (i+1,j,k)

    ixymz = ix[1:-1,0:-2,1:-1].flatten()   # index (i,j-1,k)
    ixypz = ix[1:-1,2:,1:-1].flatten()     # index (i,j+1,k)

    ixyzm = ix[1:-1,1:-1,0:-2].flatten()   # index (i,j,k-1)
    ixyzp = ix[1:-1,1:-1,2:].flatten()     # index (i,j,k+1)

    # index boundary points

    bx0=xh < h/2
    bxl=xh > L-h/2
    by0=yh < h/2
    byl=yh > L-h/2
    bz0=zh < h/2
    bzl=zh > L-h/2

    bxyz=np.logical_or.reduce((bx0,bxl,by0,byl,bz0,bzl))

    ix_BD = ix[bxyz]

    ii = np.concatenate((ixyz,ixyz,ixyz,ixyz,ixyz,ixyz,ixyz,ix_BD))     # sparse row    i
    jj = np.concatenate((ixyz,ixmyz,ixpyz,ixymz,ixypz,ixyzm,ixyzp,ix_BD)) # sparse column j

    ee  = np.ones(((k-2)**3,1))
    ee1 = np.ones((np.size(ix_BD),1))
    aa  = np.concatenate((ee*(+6)/h**3,ee*(-1)/h**3,ee*(-1)/h**3,ee*(-1)/h**3,ee*(-1)/h**3,ee*(-1)/h**3,ee*(-1)/h**3,ee1))

    Lh  = sparse.csr_matrix((aa.flatten(),(ii.flatten(),jj.flatten())),shape=(k**3,k**3))

    u0=np.zeros(k**3)
    for t in np.linspace(0,T,l):

        A=sparse.eye(k**3)+tau*Lh

        rhs = tau*np.ones(k**3).flatten()+u0
        rhs[ix_BD]=0

        u = spsolve(A, rhs)
        sns.heatmap(np.reshape(u,(k,k,k))[:,:,5],vmin=0,vmax=0.005)
        plt.show()
        u0=u
    u = np.reshape(u,(k,k,k))

    return u

'''
for slice in u:
    sns.heatmap(slice)
    plt.show()
'''
