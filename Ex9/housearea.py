#!/usr/bin/env python3
# Code written by Riccardo Parise | 412524 | Scientific Computing, M.Sc.|

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from readtria import *


def main():
    [x, y, npoint, nelement, e2p, idp, ide] = readtria("./meshes/haus")
    for k in range(nelement):
        e2 = e2p[k]
        [Fdet, Finv] = generatetransformation2D(k, e2, x, y)
    return

def generatetransformation2D(k, e2, x, y):
    zi1 = x[e2][0]
    zj1 = x[e2][1]
    zk1 = x[e2][2]
    zi2 = y[e2][0]
    zj2 = y[e2][1]
    zk2 = y[e2][2]

    A = [[zj1-zi1, zk1-zi1],[zj2-zi2, zk2-zi2]]

    print(A)
    Fdet = np.linalg.det(T)

    return


if __name__ == '__main__':
    main()
