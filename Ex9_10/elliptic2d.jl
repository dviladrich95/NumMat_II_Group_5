## FEM Julia sample code
using SparseArrays
using PyPlot

# pyplot
include("fem2d.jl") # here you can put, what you implemented

x, y, np, ne, e2, idp, ide = readtria("disc");  # read mesh from file
localtoglobal2DP1 = e2;                   # could it be this simple? why?
int = collect(.~(idp[1:np] .== 1));       # select points without Dirichlet bc
nphi = 3;

## build matrices
ii = zeros(Int64, ne, nphi, nphi); # sparse i-index
jj = zeros(Int64, ne, nphi, nphi); # sparse j-index
aa = zeros(ne, nphi, nphi); # entry of Galerkin matrix
bb = zeros(ne, nphi, nphi); # entry in mass-matrix (to build rhs)

## build global from local
for k = 1:ne             # loop over elements
    Fdet, Finv = generatetransformation2D(k, e2, x, y); # compute trafo
    
    # build local matrices (mass, stiffness, ...)
    sloc = localstiff2D(Fdet, Finv); # element stiffness matrix
    mloc = localmass2D(Fdet);       # element mass matrix
    
    # compute i,j indices of the global matrix
    dofs = localtoglobal2DP1[k,:];
   
    # compute a(i,j) values of the global matrix
    for i = 1:nphi
        for j = 1:nphi
            ii[k,i,j] = dofs[i]; # local-to-global
            jj[k,i,j] = dofs[j]; # local-to-global
            aa[k,i,j] = sloc[i,j];
            bb[k,i,j] = mloc[i,j];
        end
    end
end
# create sparse matrices
A = sparse(ii[:], jj[:], aa[:]);
M = sparse(ii[:], jj[:], bb[:]);

# build rhs and take into account Dirichlet bcs, solve, plot
rhs = M * ones(np, 1);
u   = zeros(np, 1);
u[int] = A[int, int] \ rhs[int];

# plotting
tripcolor(x[:], y[:], triangles = e2 .- 1, u[:])

