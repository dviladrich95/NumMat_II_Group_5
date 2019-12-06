#!/usr/bin/env python3
import numpy as np
import scipy.sparse
from dill.source import getsource
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------\
# Assignment 4, Exercise 1c                                               |
#                                                             submitted by|
#                                                                         |
#                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
#                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
#                   Riccardo Parise | 412524 | Scientific Computing, M.Sc.|
#        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
#                                                                         |
#                                                          in  Python 3.5 |
# ------------------------------------------------------------------------/


def main():
    # Solving for the fuctions on exercise b
    eps = [1e-1, 1e-2, 1e-3, 1e-4]
    sns.set(style="white", palette="muted", color_codes=True)

    # Setting subplots and titles
    fig, axes = plt.subplots(ncols=len(eps), nrows=3, sharex='col', sharey='row')
    row_names = ['f(x)', 'k(x)', 'u_h']
    col_names = eps
    for row in range(len(row_names)):
        axes[row, 0].axes.set_ylabel(row_names[row])
    for col in range(len(col_names)):
        axes[0, col].title.set_text("eps: " + str(col_names[col]))

    # Run the calculations for each eps
    for exp in range(len(eps)):
        def f(x): return 0 * x + 1

        def k(x): return 2 + np.tanh(1 / eps[exp] * (x - 1 / 2))

        def g(x): return 0 * x

        n = int(20 + 10 * np.sqrt(1 / eps[exp]))  # number of grid points. empirically chosen
        # function
        [Lh, xh, fh] = a06ex03getPDE(k, f, g, n)
        u_h = a06ex03solvePDE(Lh, fh)

        sns.lineplot(xh, f(xh), ax=axes[0, exp])
        sns.lineplot(xh, k(xh), ax=axes[1, exp])
        sns.lineplot(xh, u_h, ax=axes[2, exp])

    sns.despine()
    plt.show()


def a06ex03getPDE(k, f, g, N):
    print("INPUTS:")
    print("function k(x):\n", getsource(k))
    print("function f(x):\n", getsource(f))
    print("function g(x):\n", getsource(g))
    print("number of grid points N:\n", N)

    print("\nPRE-CALCULATIONS:")
    h = 1 / (N + 1)
    xh = h * np.asarray(range(1, N + 1))  # grid domain xh
    print("h: step size\n", h)
    print("xh: domain\n", xh)

    dd = scipy.sparse.diags([1, -2, 1], [-1, 0, 1], shape=[N, N])
    print("dd: standard stencil matrix d+d-\n", dd.todense())
    d0 = scipy.sparse.diags([-1 / 2, 0, 1 / 2], [-1, 0, 1], shape=[N, N])
    print("d0: standard stencil matrix d0\n", d0.todense())

    # Using the multiplication rule
    # -(k(x)*u'(x))' = f(x)
    # -k(x)*u"(x) - k'(x)*u'(x) + 0*u(x) = f(x)
    print("\nCALCULATIONS:\n")

    fh = f(xh)
    print("fh: vectors of inhomogeneity f(x):\n", fh)

    kx = k(xh)  # vector with k(x)
    mk = scipy.sparse.diags([kx], [0], shape=[N, N])
    print("mk: matrix of coefficients k(x)\n", mk.todense())

    dk = 1 / h * d0 * kx  # vector with derivative of k(x)
    dk[0] = dk[0] + k(0) * (1 / 2 / h)  # correct the value in the boundary
    dk[-1] = dk[-1] + k(1) * (1 / 2 / h)  # correct the value in the boundary
    mkdx = scipy.sparse.diags([dk], [0], shape=[N, N])
    print("mkdx: matrix of coefficients k'(x) (mkdx = 1/h*d0 * kx + BC):\n", mkdx.todense())

    Lh = 1 / (h ** 2) * -mk * dd - 1 / h * mkdx * d0
    # since the boundary condition u(0) and u(1) are 0, we do not need to deduce from the RHS fh
    print("Lh: discretization matrix (lh= 1/(h**2)*-mk*dd-1/h*mkdx*d0 + BC\n", (mk * dd).todense())
    return [Lh, xh, fh]


def a06ex03solvePDE(Lh, fh):
    u_h = scipy.sparse.linalg.spsolve(Lh, fh)
    print("u_h: solution vector of function u\n", u_h)
    return u_h


if __name__ == '__main__':
    main()