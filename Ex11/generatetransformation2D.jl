# Sample solution DP
function  generatetransformation2D(k::Int64, e2::Array{Int64,2}, x::Array{Float64,2}, y::Array{Float64,2})::Tuple{Float64,Array{Float64,2}}

    dx1 = x[e2[k,2]] - x[e2[k,1]];
    dy1 = y[e2[k,2]] - y[e2[k,1]];

    dx2 = x[e2[k,3]] - x[e2[k,1]];
    dy2 = y[e2[k,3]] - y[e2[k,1]];

# determinant on each triangle
    Fdet = dx1 * dy2 - dx2 * dy1;

# transformation jacobian on each triangle
    Finv = zeros(2, 2);
    Finv[1,1] =  dy2 / Fdet ;
    Finv[1,2] = -dx2 / Fdet ;
    Finv[2,1] = -dy1 / Fdet ;
    Finv[2,2] =  dx1 / Fdet ;
    return Fdet, Finv
end