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

def generatetransformation2D(k, e2, x, y):
    i, j, k = e2[k, :]
    zi1, zj1, zk1 = x[i], x[j], x[k]
    zi2, zj2, zk2 = y[i], y[j], y[k]

    F_l = np.array([[zj1 - zi1, zk1 - zi1], [zj2 - zi2, zk2 - zi2]])

    Fdet = np.linalg.det(F_l)

    Finv = np.linalg.inv(F_l)

    return Fdet, Finv


# matrices from previous assignment
def localmass2D(Fdet):
    mref = np.array([[1 / 12, 1 / 24, 1 / 24], [1 / 24, 1 / 12, 1 / 24], [1 / 24, 1 / 24, 1 / 12]])
    mloc = mref * abs(Fdet)
    return mloc


def localstiff2D(Fdet, Finv):
    B_l = Finv.dot(Finv.transpose())
    grad_phi = np.array([[-1, 1, 0], [-1, 0, 1]])
    sloc = abs(Fdet)/2 * np.linalg.multi_dot([grad_phi.transpose(), B_l, grad_phi])
    return sloc



if __name__ == "__main__":
    # FEM Python sample code
    x, y, npo, ne, e2, idp, ide = readtria('./meshes/disc')  # read mesh from file
    localtoglobal2DP1 = e2  # could it be this simple? why?
    # select points without Dirichlet bc
    it = np.logical_not(idp == 1)
    nphi = 3

    # build matrices
    ii = np.zeros((ne, nphi ** 2))  # sparse i-index
    jj = np.zeros((ne, nphi ** 2))  # sparse j-index
    aa = np.zeros((ne, nphi ** 2))  # entry of Galerkin matrix
    bb = np.zeros((ne, nphi ** 2))  # entry in mass-matrix (to build rhs)

    for k in np.arange(0, ne):
        Fdet, Finv = generatetransformation2D(k, e2, x, y)  # compute trafo

        # build local matrices (mass, stiffness, ...)
        sloc = localstiff2D(Fdet, Finv)  # element stiffness matrix
        mloc = localmass2D(Fdet)  # element mass matrix

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
    rhs = M * np.ones(npo)
    u = np.zeros(npo)
    u[it] = spsolve(A[np.ix_(it, it)], rhs[it])

    # plotting
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # , linewidth=0.2, antialiased=True)
    u_exact = -1 / 8 * (x ** 2 + y ** 2) + 0.5
    ax.plot_trisurf(x, y, u, triangles=e2, alpha=0.6)
    ax.plot_trisurf(x, y, u_exact, alpha=0.3)
    plt.show()
    print("Error is %f", sum(u - u_exact) / len(u))
