% parameters
clear all
close all

N = 1022;
L = 1;
k = N + 2;    % define k = N+2 for ease of typing
h = L/(k-1);  % lattice spacing

% array construction
[xh,yh]=meshgrid(linspace(0,L,k)); % create mesh using MATLABs meshgrid
xh = xh';x = xh(:);
yh = yh';y = yh(:);

ix = reshape(1:(k^2),[k k]);    % create indices in lexic. ordering

ixy               = ix((2:k-1)  ,(2:k-1)  );ixy =ixy(:);  % index (i,j)
ixmy(1:k-2,1:k-2) = ix((2:k-1)-1,(2:k-1)  );ixmy=ixmy(:); % index (i-1,j)
ixpy(1:k-2,1:k-2) = ix((2:k-1)+1,(2:k-1)  );ixpy=ixpy(:); % index (i+1,j)
ixym(1:k-2,1:k-2) = ix((2:k-1)  ,(2:k-1)-1);ixym=ixym(:); % index (i,j-1)
ixyp(1:k-2,1:k-2) = ix((2:k-1)  ,(2:k-1)+1);ixyp=ixyp(:); % index (i,j+1)

ix_BD = ix( xh < h/2 | yh < h/2 | xh > L-h/2 | yh > L-h/2 ); % index boundary points
ix_BD = ix_BD(:);

% Q: How to choose aa, so that Lh becomes the 5-point stencil
% approximation of the Laplacian using the code below?
%{
% Some visual aid in order to check the functionality of the arrays
% for small N.
%
plot(x,y,'k.',x(ix_BD),y(ix_BD),'bo',x(ixy),y(ixy),'rx');
for i=1:length(x)
    text(x(i),y(i)+h/8,sprintf('%i',ix(i)))
end
xlim([-0.1 1.1])
ylim([-0.1 1.1])
%}

%% SOLUTION (without function)
%{-
ii = [ixy;ixy;ixy;ixy;ixy;ix_BD];     % sparse row    i
jj = [ixy;ixmy;ixpy;ixym;ixyp;ix_BD]; % sparse column j
aa = [ones((k-2)^2,1) *(+4)/h^2;... % u(x,y)
    ones((k-2)^2,1) *(-1)/h^2;... % u(x-1,y)
    ones((k-2)^2,1) *(-1)/h^2;... % u(x+1,y)
    ones((k-2)^2,1) *(-1)/h^2;... % u(x,y-1)
    ones((k-2)^2,1) *(-1)/h^2;... % u(x,y+1)
    ones(length(ix_BD),1)];

Lh = sparse(ii,jj,aa);

rhs = sin(pi*x).*sin(pi*y);
u   = Lh\rhs;
u   = reshape(u,[N+2 N+2]);

surf(xh,yh,u, 'EdgeAlpha', 0.1);
title(sprintf('dofs %i, error %e',(N+2)^2,max(max(abs(u-sin(pi*xh).*sin(pi*yh)/(2*pi^2))))));
%}






