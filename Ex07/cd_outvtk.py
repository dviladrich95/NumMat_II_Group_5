from a07ex02getsola import a07ex02getsola
from a07ex02b import a07ex02b

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


R   = 10

N=30

L=1.0

Nx=N
Ny=N
Nz=N

ua = a07ex02getsola(Nx,Ny,Nz,R)
ub = a07ex02b(N)

udiff=ua-ub

outvtk_structured3d('a07ex02vtk_a',L,N,ua.flatten(order="F"))
outvtk_structured3d('a07ex02vtk_b',L,N,ub.flatten(order="F"))
outvtk_structured3d('a07ex02vtk_diff',L,N,udiff.flatten(order="F"))




