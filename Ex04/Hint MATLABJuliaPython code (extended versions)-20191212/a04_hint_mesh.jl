# parameters
using SparseArrays
using PyPlot

N = 128;
L = 1;
k = N + 2;    # define k = N+2 for ease of typing
h = L / (k - 1);  # lattice spacing

function meshgrid(x::LinRange{Float64})
    N = length(x);
    xh = zeros(N, N);
    yh = zeros(N, N);
    for i = 1:N
        for j = 1:N
            xh[i,j] = x[i];
            yh[i,j] = x[j];
        end
    end
    return xh, yh
end

# array construction
xh, yh = meshgrid(LinRange(0, L, k)); # create mesh using MATLABs meshgrid

x = xh[:];
y = yh[:];

ix = reshape(Vector(1:k^2), (k, k));

ixy  = ix[(2:k - 1),(2:k - 1)]; # index (i,j)
ixmy = ix[(1:k - 2),(2:k - 1)]; # index (i-1,j)
ixpy = ix[(3:k    ),(2:k - 1)]; # index (i+1,j)
ixym = ix[(2:k - 1),(1:k - 2)]; # index (i,j-1)
ixyp = ix[(2:k - 1),(3:k    )]; # index (i,j+1)

ix_BD = @. ix[ (xh < h / 2) | (yh < h / 2) | (xh > L - h / 2) | (yh > L - h / 2)]; # index boundary points

ixy  = ixy[:];
ixmy = ixmy[:];
ixpy = ixpy[:];
ixym = ixym[:];
ixyp = ixyp[:];


# Q: How to choose aa, so that Lh becomes the 5-point stencil
# approximation of the Laplacian using the code below?

# Some visual aid in order to check the functionality of the arrays
# for small N.
# scatter(x, y)
# scatter(x[ix_BD], y[ix_BD], color = :red)

ii = [ixy;ixy;ixy;ixy;ixy;ix_BD];     # sparse row    i
jj = [ixy;ixmy;ixpy;ixym;ixyp;ix_BD]; # sparse column j
aa = [ones((k - 2)^2) * (+4) / h^2; # u(x,y)
    ones((k - 2)^2) * (-1) / h^2; # u(x-1,y)
    ones((k - 2)^2) * (-1) / h^2; # u(x+1,y)
    ones((k - 2)^2) * (-1) / h^2; # u(x,y-1)
    ones((k - 2)^2) * (-1) / h^2; # u(x,y+1)
    ones(length(ix_BD))];

Lh = sparse(ii, jj, aa);

rhs = @. sin(pi * x) .* sin(pi * y);
u   = Lh \ rhs;
u   = reshape(u, (k, k));

surf(xh, yh, u, cmap = ColorMap("jet"))
println("dofs:", (N + 2)^2);
println("error:", maximum(abs.(u - sin.(pi * xh) .* sin.(pi * yh) / (2 * pi^2))));