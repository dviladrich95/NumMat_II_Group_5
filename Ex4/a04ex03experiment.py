import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt

def analytical_u(eps, xh):
    # Analitical solution for the equation using eps as variable.

    if eps != 0:
        ua = list(x - (np.exp(-(1 - x) / eps) - np.exp(-1 / eps)) / (1 - np.exp(-1 / eps)) for x in xh)
    else:
        # Solution for lim eps -> 0
        ua = xh
    return ua

def a04ex03solve(eps,xh,flag):
    # Create required matrix for the equation using epsm and solve using linealg library.
    # The first derivative can be chosen between D-, D+ and D0 using the flag +, -, 0.

    h = np.diff(xh)         # step size
    N = len(h) - 1          # number of equation (not included the BC)
    fh = np.ones(len(h)-1)  # RHS of equation


    f_dd_m1 = lambda i: 2 / (h[i] * (h[i] + h[i + 1]))
    f_dd_0 = lambda i: -2 / (h[i] * h[i + 1])
    f_dd_p1 = lambda i: 2 / (h[i + 1] * (h[i] + h[i + 1]))

    f_d_0_1 = lambda i: -1 / (h[i] + h[i + 1])
    f_d_0_2 = lambda i: 1 / (h[i] + h[i + 1])

    f_d_p_1 = lambda i: -1 / h[i + 1]
    f_d_p_2 = lambda i: 1 / h[i + 1]

    f_d_m_1 = lambda i: -1 / h[i]
    f_d_m_2 = lambda i: 1 / h[i]

    # Second derivative matrix
    dd_m1 = [f_dd_m1(i) for i in range(1, N)]
    dd_0 = [f_dd_0(i) for i in range(N)]
    dd_p1 = [f_dd_p1(i) for i in range(N - 1)]
    dd = np.diag(dd_m1, k=-1) + np.diag(dd_0, k=0) + np.diag(dd_p1, k=1)

    # First derivartive matrix
    if flag == '-':
        d_1 = [f_d_m_1(i) for i in range(1, N)]
        d_2 = [f_d_m_2(i) for i in range(N)]

        d = np.diag(d_1, k=-1) + np.diag(d_2, k=0)
    elif flag == '+':
        d_1 = [f_d_p_1(i) for i in range(N)]
        d_2 = [f_d_p_2(i) for i in range(N - 1)]
        d = np.diag(d_1, k=0) + np.diag(d_2, k=1)
    elif flag == '0':
        d_1 = [f_d_0_1(i) for i in range(1, N)]
        d_2 = [f_d_0_2(i) for i in range(N - 1)]
        d = np.diag(d_1, k=-1) + np.diag(d_2, k=1)


    # Set up the equation in the matrix form
    l_h = sps.csr_matrix(-eps * dd + d)

    # Solve
    uh = sps.linalg.spsolve(l_h, fh)

    # Add the boundary conditions
    uh = np.concatenate(([0], uh, [0]))
    return uh


def a04ex03error(eps,xh,uh):
    # Calculate the analytical solution for the grip points xh using the eps, and compare the error with the
    # numerical solution uh

    uex = analytical_u(eps,xh)  # Solve the equation for the grip points
    err = max(abs(uex - uh))    # Get norm inf of the difference
    return err, uex


def a04ex03shishkin(N, sigma):
    # Generates the Shishkin mesh
    xh_1 = list(i*(1 - sigma) / N for i in range(N))
    xh_2 = list((1 - sigma) + (i - N ) * sigma / N for i in range(N+1, 2*N + 1))
    xh = xh_1 + xh_2
    return xh


## Program Start

# User defined Variables
eps = 0.001
exps = ['i', 'ii', 'iii', 'iv']   # Experiment title
Ns = [5, 50, 500, 5000]           # Mesh size
exps_mesh = ['u', 'u', 'u', 's']  # Mesh type:  'u' = uniform, s = shishkin
exps_flag = ['+', '0', '-', '0']  # First derivative flag:  '+' = D^+, '-' = D^-, '0' = D^0,


for i in range(len(exps)):
    sigma = 4 * eps * np.log(2 * Ns[i])

    # Initialize plot
    fig = plt.figure()
    fig.suptitle('Experiment ' + exps[i] + ') - Mesh type: ' + exps_mesh[i] + '   Flag: ' + exps_flag[i])
    fig.subplots_adjust(hspace=0.5, wspace=0.5)

    for i in range(4):
        # Mesh Generation depending on the flag exps_mesh
        xa = np.linspace(0, 1, 1000)                   # Mesh for smooth plot
        if exps_mesh[i] == 'u':
            xh = np.linspace(0, 1, num=2 * Ns[i] + 1)  # Uniform
        elif exps_mesh[i] == 's':
            xh = a04ex03shishkin(Ns[i], sigma)         # Shishkin mesh

        # Calculations
        ua = analytical_u(eps, xa)                # Calculate the analytical solution for the fine mesh to plot
        uh = a04ex03solve(eps, xh, exps_flag[i])  # Calculate the numerical solution for the mesh xh
        [err, uex] = a04ex03error(eps, xh, uh)    # Calculate the error for the solution uh

        # Plots
        ax = fig.add_subplot(2, 2, i+1)
        ax = sns.lineplot(xa, ua, color='black')
        ax = sns.scatterplot(xh, uh, linewidth=0.1, edgecolor='black', color=(0.333, 0.659, 0.408))
        ax.set(xlabel='x', ylabel='u')
        ax.title.set_text('N = %d, Err = %3.1e' % (Ns[i], err))
plt.show()