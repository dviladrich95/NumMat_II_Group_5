import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt


def generateelements1D(x):
    npo = x.size                            # number of points
    ne = npo-1                              # number of elements
    a = np.arange(0, npo-1).reshape((-1, 1))
    b = np.arange(1, npo).reshape((-1, 1))
    e1 = np.concatenate((a, b), axis=1)      # element specificiation
    return ne, npo, e1


def generatetransformation1D(k, e1, x):
    dx = x[e1[k, 1]]-x[e1[k, 0]]            # length of k-th interval
    Fdet = dx                             # determinant
    Finv = 1/dx                             # inverse transformation
    return Fdet, Finv


def localmass1D(Fdet):
    return np.array(((1./3, 1./6), (1./6, 1./3)))*Fdet


def localstiff1D(Fdet, Finv):
    G = np.array(Finv)
    sloc = np.zeros((2, 2))
    sloc[0, 0] = G.T.dot(G.dot(Fdet))
    sloc[0, 1] = - G.T.dot(G.dot(Fdet))
    sloc[1, 0] = sloc[0, 1]
    sloc[1, 1] = sloc[0, 0]
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
    A = sparse.csr_matrix((aa.flatten(order="F"), (ii.flatten(
        order="F"), jj.flatten(order="F"))), shape=(N, N))
    M = sparse.csr_matrix((bb.flatten(order="F"), (ii.flatten(
        order="F"), jj.flatten(order="F"))), shape=(N, N))

    # build rhs and take into account Dirichlet bcs, solve, plot
    rhs = M.dot(np.ones(N))
    u = np.zeros(N)
    u[1:-1] = spsolve(A[1:-1, 1:-1], rhs[1:-1])

    plt.figure()
    plt.plot(x, u, x, x*(1-x)/2, 'r.')
    plt.show()
