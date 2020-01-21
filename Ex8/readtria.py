#!/usr/bin/python3
# Note: Code probably needs to be changed if extra marker are introduced
import numpy as np
import matplotlib.pyplot as plt

def readtria(fname):
    # input: filename
    # output: x,y,npoint,nelement,e2p,idp,ide
    print("read mesh file: " + fname)  
    
    # Node file
    npoint = ndim = npattr = nbound = None    
    fname_node = fname + ".node"
    with open(fname_node, "r") as fnode:    
        npoint, ndim, npattr, nbound = [int(x) for x in fnode.readline().split()]         

    U = np.genfromtxt(fname_node, skip_header=1)
    x = U[:,1]               # dimension
    y = U[:,2]               # number of point attributes
    idp = U[:,3].astype(int) # number of boundary markers

    
    
    # Element file
    nelement = nphi = neattr = None    
    fname_ele = fname + ".ele"
    with open(fname_ele, "r") as fele:
        nelement, nphi, neattr = [int(x) for x in fele.readline().split()]         

    U = np.genfromtxt(fname_ele, skip_header=1, dtype=np.int64)
    e2p = U[:, 1:nphi+1] -1
    ide = U[:, nphi+1:]
    
    return x,y,npoint,nelement,e2p,idp,ide


x,y,npoint,nelement,e2p,idp,ide = readtria("box.2")
plt.scatter(x,y)
for tri in e2p:
    plt.plot([x[tri[0]], x[tri[1]]], [y[tri[0]], y[tri[1]]])
    plt.plot([x[tri[1]], x[tri[2]]], [y[tri[1]], y[tri[2]]])
    plt.plot([x[tri[0]], x[tri[2]]], [y[tri[0]], y[tri[2]]])
plt.show()
print("Number of triangular elements: ", int(len(e2p)/3))
print("Number of nodes: ", len(x))