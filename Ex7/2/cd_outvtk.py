import numpy as np
from a07ex02getsola import a07ex02getsola
from a07ex02b import func4b

def outvtk_structured3d(fname, L, N, uh):
    h = L / (N + 1)
    uh = uh.flatten(order="F")

    with open(fname + ".vtk", "w") as f:
        f.write('# vtk DataFile Version 2.0\n')
        f.write('Num2Ing Structured 3D data\n')
        f.write('ASCII\n')
        f.write('DATASET STRUCTURED_POINTS\n')
        f.write('DIMENSIONS {N} {N} {N}\n'.format(N=N + 2))
        f.write('ORIGIN 0.0 0.0 0.0\n')
        f.write('SPACING {h} {h} {h}\n'.format(h=h))
        f.write('POINT_DATA {0}\n'.format((N + 2) ** 3))
        f.write('SCALARS solution float 1\n')
        f.write('LOOKUP_TABLE default\n')
        # f.write("\n".join([str(x) for x in uh]) + "\n")
        f.write("\n".join(["{:.6f}".format(x) for x in uh]) + "\n")


R   = 1

L=1

N=30

Nx=N
Ny=N
Nz=N

x=np.linspace(0,1,Nx+2)
y=np.linspace(0,1,Ny+2)
z=np.linspace(0,1,Nz+2)

ua = a07ex02getsola(x,y,z,R)
ub = func4b()

udiff=ua-ub
outvtk_structured3d('a07ex02vtk_a',L,N,ua.flatten(order="F"))
outvtk_structured3d('a07ex02vtk_b',L,N,ub.flatten(order="F"))
outvtk_structured3d('a07ex02vtk_diff',L,N,udiff.flatten(order="F"))




