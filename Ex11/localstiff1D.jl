# Sample solution DP
function localstiff1D(Fdet::Float64, Finv::Float64)::Array{Float64,2}
    G = Finv;
    sloc = zeros(2, 2);
    sloc[1,1] =   G' * G * Fdet;
    sloc[1,2] = - G' * G * Fdet;
    sloc[2,1] = sloc[1,2];
    sloc[2,2] = sloc[1,1];
    return sloc;
end