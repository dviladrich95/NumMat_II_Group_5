% Sample solution DP
function [Fdet,Finv]=generatetransformation1D(k,e1,x)

dx   = x(e1(k,2))-x(e1(k,1)); % length of k-th interval
Fdet =   dx; % determinant
Finv = 1/dx; % inverse transformation
