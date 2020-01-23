% Local P1 stiffness matrix
function S=localstiff2D(Fdet,Finv)
gradphi = [-1 -1;1 0;0 1];
dphi    = gradphi*Finv;
S=1/2*(dphi(:,1)*dphi(:,1)'+dphi(:,2)*dphi(:,2)')*Fdet;
