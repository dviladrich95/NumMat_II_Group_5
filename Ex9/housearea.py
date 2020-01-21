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
    T = np.zeros([3, 3])
    T[2, :] = y[e2]
    Fdet = np.linalg.det(T)

    return


if __name__ == '__main__':
    main()w
