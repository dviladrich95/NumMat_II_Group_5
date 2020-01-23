# Local P1 mass matrix
function localmass2D(Fdet::Float64)::Array{Float64,2}
    mloc = Fdet * [1 1 / 2 1 / 2;1 / 2 1 1 / 2;1 / 2 1 / 2 1] / 12;
    return mloc;
end