import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from readtria import readtria

# These you need to implement:
# def generatetransformation2D(k,e2,x,y):
# def localstiff2D(Fdet,Finv):
# def localmass2D(Fdet):


if __name__ == "__main__":
    # FEM Python sample code
    x, y, npo, ne, e2, idp, ide = readtria('box')  # read mesh from file
    localtoglobal2DP1 = e2                  # could it be this simple? why?
    # select points without Dirichlet bc
    it = np.logical_not(idp == 1)
    nphi = 3

    # build matrices
    ii = np.zeros((ne, nphi**2))  # sparse i-index
    jj = np.zeros((ne, nphi**2))  # sparse j-index
    aa = np.zeros((ne, nphi**2))  # entry of Galerkin matrix
    bb = np.zeros((ne, nphi**2))  # entry in mass-matrix (to build rhs)

    for k in np.arange(0, ne):
        Fdet, Finv = generatetransformation2D(k, e2, x, y)  # compute trafo

        # build local matrices (mass, stiffness, ...)
        sloc = localstiff2D(Fdet, Finv)  # element stiffness matrix
        mloc = localmass2D(Fdet)       # element mass matrix

        # compute i,j indices of the global matrix
        dofs = localtoglobal2DP1[k, :]
        ii[k, :] = dofs[[0, 1, 2, 0, 1, 2, 0, 1, 2]]  # local-to-global
        jj[k, :] = dofs[[0, 0, 0, 1, 1, 1, 2, 2, 2]]  # local-to-global

        # compute a(i,j) values of the global matrix
        aa[k, :] = sloc.flatten(order="F")
        bb[k, :] = mloc.flatten(order="F")

    # create sparse matrices
    A = sparse.csr_matrix((aa.flatten(order="F"), (ii.flatten(
        order="F"), jj.flatten(order="F"))), shape=(npo, npo))
    M = sparse.csr_matrix((bb.flatten(order="F"), (ii.flatten(
        order="F"), jj.flatten(order="F"))), shape=(npo, npo))

    # build rhs and take into account Dirichlet bcs, solve, plot
    rhs = M*np.ones(npo)
    u = np.zeros(npo)
    u[it] = spsolve(A[np.ix_(it, it)], rhs[it])

    # plotting
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # , linewidth=0.2, antialiased=True)
    ax.plot_trisurf(x, y, u, triangles=e2)
    plt.show()
