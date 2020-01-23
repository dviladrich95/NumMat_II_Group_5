# Local P1 stiffness matrix
function localstiff2D(Fdet::Float64, Finv::Array{Float64,2})::Array{Float64,2}
    gradphi = [-1 -1;1 0;0 1];
    dphi    = gradphi * Finv;
    sloc = 1 / 2 * (dphi[:,1] * dphi[:,1]' + dphi[:,2] * dphi[:,2]') * Fdet;
    return sloc;
end
