# Sample solution DP
function generatetransformation1D(k::Int64, e1::Array{Int64,2}, x::LinRange{Float64})::Tuple{Float64,Float64}
    dx   = x[e1[k,2]] - x[e1[k,1]]; # length of k-th interval
    Fdet =   dx; # determinant
    Finv = 1 / dx; # inverse transformation
    return Fdet, Finv;
end