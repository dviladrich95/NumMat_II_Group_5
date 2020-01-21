#!/usr/bin/env python3
# Code written by Riccardo Parise | 412524 | Scientific Computing, M.Sc.|

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from readtria import *

def main():
    # Exercise b, check generatetransformation2D by calculating the area of the house mesh
    [x, y, npoint, nelement, e2p, idp, ide] = readtria("./meshes/haus")
    # initialize area
    area = 0
    # add each element
    for k in range(nelement):
        [Fdet, Finv] = generatetransformation2D(k, e2p, x, y)
        area = area+Fdet/2  # have to divide by two. area of reference element is 1/2
    print(area)
    return

def generatetransformation2D(k, e2p, x, y):
    e2 = e2p[k]
    zi1 = x[e2][0]
    zj1 = x[e2][1]
    zk1 = x[e2][2]
    zi2 = y[e2][0]
    zj2 = y[e2][1]
    zk2 = y[e2][2]

    Fdet = np.linalg.det([[zj1-zi1, zk1-zi1],[zj2-zi2, zk2-zi2]])
    Finv = 1/Fdet

    return Fdet, Finv

def localmass2D(Fdet):
    mloc = 1
    return mloc


def localstiff2D(Fdet,Finv):
    sloc = 1
    return sloc


if __name__ == '__main__':
    main()
