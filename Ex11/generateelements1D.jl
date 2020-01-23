# Sample solution DP
function generateelements1D(x::LinRange{Float64})::Tuple{Int64,Int64,Array{Int64,2}}
    np = length(x);         # number of points
    ne = np - 1;              # number of elements
    e1 = [collect(1:(np - 1)) collect(2:np)];  # element specificiation
    return ne, np, e1;
end