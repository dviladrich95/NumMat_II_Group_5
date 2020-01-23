## FEM Julia sample code
using SparseArrays
using Plots
include("fem1d.jl") # here you can put, what you implemented

N    =  32;                         # 1D: number of points
nphi =   2;                         # P1: 2 local basis functions
x = LinRange(0, 1, N);              # array of points
ne, np, e1 = generateelements1D(x); # generate mesh
localtoglobal1DP1 = e1;             # could it be this simple? why?

## build matrices
ii = zeros(Int64, ne, nphi, nphi); # sparse i-index
jj = zeros(Int64, ne, nphi, nphi); # sparse j-index
aa = zeros(ne, nphi, nphi); # entry of Galerkin matrix
bb = zeros(ne, nphi, nphi); # entry in mass-matrix (to build rhs)

## build global from local
for k = 1:ne             # loop over elements
    Fdet, Finv = generatetransformation1D(k, e1, x); # compute trafo
    
    # build local matrices (mass, stiffness, ...)
    sloc = localstiff1D(Fdet, Finv); # element stiffness matrix
    mloc = localmass1D(Fdet);       # element mass matrix
    
    # compute i,j indices of the global matrix
    dofs = localtoglobal1DP1[k,:];
   
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
rhs = M * ones(N, 1);
u   = zeros(N, 1);
u[2:end - 1] = A[2:end - 1,2:end - 1] \ rhs[2:end - 1];

xh = collect(x);
plot(xh, u)
scatter!(xh, xh .* (1 .- xh) ./ 2)
