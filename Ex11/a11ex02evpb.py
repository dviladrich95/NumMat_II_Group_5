import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve, eigs, eigsh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from readtria import readtria

def generatetransformation2D(k, e2, x, y):
    dx1 = x[e2[k, 1]]-x[e2[k, 0]]
    dy1 = y[e2[k, 1]]-y[e2[k, 0]]

    dx2 = x[e2[k, 2]]-x[e2[k, 0]]
    dy2 = y[e2[k, 2]]-y[e2[k, 0]]

    # determinant on each triangle
    Fdet = dx1*dy2 - dx2*dy1

    # transformation jacobian on each triangle
    Finv = np.zeros((2, 2))
    Finv[0, 0] = dy2 / Fdet
    Finv[0, 1] = -dx2 / Fdet
    Finv[1, 0] = -dy1 / Fdet
    Finv[1, 1] = dx1 / Fdet

    return Fdet, Finv


def localstiff2D(Fdet, Finv):
    gradphi = np.array([[-1, -1], [1, 0], [0, 1]])
    dphi = gradphi.dot(Finv)
    S = 1/2*(dphi[:, 0].reshape(-1, 1).dot(dphi[:, 0].reshape(1, -1)) +
             dphi[:, 1].reshape(-1, 1).dot(dphi[:, 1].reshape(1, -1)))*Fdet
    return S

# Local P1 mass matrix


def localmass2D(Fdet):
    return Fdet*np.array([[1, 1/2, 1/2], [1/2, 1, 1/2], [1/2, 1/2, 1]])/12


if __name__ == "__main__":
    # FEM Python sample code
    x, y, npo, ne, e2, idp, ide = readtria('mdisc')  # read mesh from file
    # select points without Dirichlet bc
    it = np.logical_not(idp == 1)

    # fake mesh
    # x = np.array([0, 1, 1, 0, 0.5])
    # y = np.array([0, 0, 1, 1, 0.5])
    # e2 = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
    # # select points without Dirichlet bc
    # it = np.array([False, False, False, False, True])
    # ne = 4
    # npo = 5
    # idp = np.array([1, 1, 1, 1, 0])
    # ide = np.array([[1], [1], [1], [1], [1]])


    localtoglobal2DP1 = e2                  # could it be this simple? why?
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
    w, v = eigs(A[np.ix_(it, it)], M=M[np.ix_(it, it)], sigma=1)

    # plotting
    fig, ax = plt.subplots(1, 6, figsize=(18, 4))
    fig.suptitle('Eigenvalues')
    for i in range(6):
        u[it] = np.real(v[:, i])
        ax[i].tricontourf(x, y, u, 100, cmap="inferno")
        ax[i].title.set_text("%3.1f + %3.1fj" % (np.real(w[i]), np.imag(w[i])))
        ax[i].axis('off')
    plt.show()
