import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# These you need to implement:
# def generateelements1D(x):
# def generatetransformation1D(k, e1, x):
# def localmass1D(Fdet):
# def localstiff1D(Fdet, Finv):

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
