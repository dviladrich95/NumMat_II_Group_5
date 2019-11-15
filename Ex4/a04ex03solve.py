import numpy as np
import scipy.sparse as sps
import seaborn as sns
import matplotlib.pyplot as plt


def analytical_u(xh, eps):
    return (x - (np.exp(-(1 - x) / eps) - np.exp(-1 / eps)) / (1 - np.exp(-1 / eps)) for x in xh)

def numerical_u(xh, eps, flag):
    h = np.diff(xh)
    N = len(h) - 1
    f_h = np.ones(len(h)-1)

    dd_m1 = [2 / (h[i] * (h[i] + h[i + 1])) for i in range(N - 1)]
    dd_0 = [-2 / (h[i] * h[i + 1]) for i in range(N)]
    dd_p1 = [2 / (h[i + 1] * (h[i] + h[i + 1])) for i in range(N - 1)]
    dd = np.diag(dd_m1, k=-1) + np.diag(dd_0, k=0) + np.diag(dd_p1, k=1)

    if flag == '-':
        d_1 = [-1 / h[i] for i in range(N - 1)]
        d_2 = [1 / h[i] for i in range(N)]
        d = np.diag(d_1, k=-1) + np.diag(d_2, k=0)
    elif flag == '+':
        d_1 = [-1 / h[i + 1] for i in range(N)]
        d_2 = [1 / h[i + 1] for i in range(N - 1)]
        d = np.diag(d_1, k=0) + np.diag(d_2, k=1)
    elif flag == '0':
        d_1 = [-1 / (h[i] + h[i + 1]) for i in range(N - 1)]
        d_2 = [1 / (h[i] + h[i + 1]) for i in range(N - 1)]
        d = np.diag(d_1, k=-1) + np.diag(d_2, k=1)

    l_h = -eps * dd + d

    u_h = sps.linalg.spsolve(l_h, f_h)
    u_h = np.concatenate(([0], u_h, [0]))
    return u_h


# Setting mesh size
N = 10
eps = 0.1

# Making the non uniform mesh for numerical solver
xh = list(x for x in np.linspace(0, 1, num=N + 2) + np.random.uniform(-1 / (N + 2), 1 / (N + 2), N + 2))
# Making sure the B.C. are right despite the induced noise
xh[0] = 0
xh[-1] = 1
# Making finer mesh for analitical solver
xa = np.linspace(0, 1, 100)

# Get the analitical solution
u_a = analytical_u(xa, eps)

ax1 = sns.lineplot(xa, u_a)
ax1.set(xlabel='x', ylabel='u_a')

# Get the numerical solution
u_h = numerical_u(xh, eps, '-')

ax2 = sns.scatterplot(xh, u_h)
ax2.set(xlabel='x', ylabel='u_h')
plt.show()

#
# ax = sns.scatterplot(xh, analytical_u(x))
# for noise in np.linspace(0, 0.0, num=2):
#     xh = np.linspace(0, 1, num=N + 2) ** 2 + noise / 2 * np.random.uniform(-1 / (N + 2), 1 / (N + 2), N + 2)
#     xh[0] = 0
#     xh[-1] = 1
#
#     f = analytical_f(x)[1:-1]
#     flag = '-'
#     u_h = a04ex02solve(x, f, consts, flag)
#     sns.scatterplot(x, u_h)
# ax.set(xlabel='x', ylabel='u_h')
# plt.show()
