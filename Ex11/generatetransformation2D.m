% Sample solution DP
function [Fdet,Finv] = generatetransformation2D(k,e2,x,y)

dx1=x(e2(k,2))-x(e2(k,1));
dy1=y(e2(k,2))-y(e2(k,1));

dx2=x(e2(k,3))-x(e2(k,1));
dy2=y(e2(k,3))-y(e2(k,1));

% determinant on each triangle
Fdet = dx1*dy2 - dx2*dy1;

% transformation jacobian on each triangle
Finv(1,1) =  dy2 / Fdet ;
Finv(1,2) = -dx2 / Fdet ;
Finv(2,1) = -dy1 / Fdet ;
Finv(2,2) =  dx1 / Fdet ;