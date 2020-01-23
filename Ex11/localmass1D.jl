# Sample solution DP
function localmass1D(Fdet::Float64)::Array{Float64,2}
    mloc = [1 / 3 1 / 6;1 / 6 1 / 3] * Fdet;
    return mloc;
end