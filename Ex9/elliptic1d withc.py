import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt


def generateelements1D(x):
    num_e=len(x)-1
    num_p=len(x)

    i1_list=np.arange(num_e)
    i2_list=np.arange(1,num_e+1)

    e1=np.stack([i1_list,i2_list])
    e1=e1.transpose()
    return num_e, num_p, e1

def generatetransformation1D(k,e1,x):
    x_bar=x#np.concatenate([[0],x,[1]])
    i1,i2=e1[k,:]
    a,b=x_bar[i1],x_bar[i2]
    Fdet = a-b
    Finv=1/(b-a)
    return Fdet, Finv

#matrices from previous assignment
def localmass1D(Fdet):
    mloc=1/6*np.array([[2,1],[1,2]])*abs(Fdet)
    return mloc
def localstiff1D(Fdet,Finv):
    sloc=np.array([[1,-1],[-1,1]])*Finv**2*abs(Fdet)
    return sloc


if __name__ == "__main__":
    # FEM MATLAB sample code
    N = 32                         # 1D: number of points
    nphi = 2                         # P1: 2 local basis functions
    x = np.linspace(0, 1, N)             # array of points
    ne, npo, e1 = generateelements1D(x)  # generate mesh
    localtoglobal1DP1 = e1             # could it be this simple? why?
    # build matrices
    ii = np.zeros((ne, nphi**2))  # sparse i-index
    jj = np.zeros((ne, nphi**2))  # sparse j-index
    aa = np.zeros((ne, nphi**2))  # entry of Galerkin matrix
    bb = np.zeros((ne, nphi**2))  # entry in mass-matrix (to build rhs)
    for k in np.arange(0, ne):
        Fdet, Finv = generatetransformation1D(k, e1, x)  # compute trafo

        # build local matrices (mass, stiffness, ...)
        sloc = localstiff1D(Fdet, Finv)  # element stiffness matrix
        mloc = localmass1D(Fdet)       # element mass matrix

        # compute i,j indices of the global matrix
        dofs = localtoglobal1DP1[k, :]
        ii[k, :] = [dofs[0], dofs[1], dofs[0], dofs[1]]  # local-to-global
        jj[k, :] = [dofs[0], dofs[0], dofs[1], dofs[1]]  # local-to-global

        # compute a(i,j) values of the global matrix
        aa[k, :] = sloc.flatten(order="F")
        bb[k, :] = mloc.flatten(order="F")

    # create sparse matrices
    A = sparse.csr_matrix((aa.flatten(order="F"), (ii.flatten(order="F"), jj.flatten(order="F"))), shape=(N, N))
    M = sparse.csr_matrix((bb.flatten(order="F"), (ii.flatten(order="F"), jj.flatten(order="F"))), shape=(N, N))
    A=A+M
    # build rhs and take into account Dirichlet bcs, solve, plot
    rhs = M.dot(np.ones(N))
    u = np.zeros(N)
    u[1:-1] = spsolve(A[1:-1, 1:-1], rhs[1:-1])

    plt.figure()
    f_with_c=np.exp(-x)*(1-np.exp(x))*(np.exp(x)-np.e)/(1+np.e)
    plt.plot(x, u, x,f_with_c, 'r.')
    #plt.plot(x, u, x, x*(1-x)/2, 'r.')
    plt.show()
