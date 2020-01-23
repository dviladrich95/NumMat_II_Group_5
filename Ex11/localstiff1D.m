% Sample solution DP
function sloc=localstiff1D(Fdet,Finv)
G = Finv;
sloc(1,1) =   G' * G * Fdet;
sloc(1,2) = - G' * G * Fdet;
sloc(2,1) = sloc(1,2);
sloc(2,2) = sloc(1,1);
end