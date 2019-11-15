function [Lh, Xh, Yh, Idx_Bd] = a04ex01_Lh5(L, N)

h = L / (N + 1);  % lattice spacing
[Xh, Yh]=meshgrid(linspace(0, L, N + 2)); % create mesh using MATLABs meshgrid

Xh = Xh';
X  = Xh(:);

Yh = Yh';
Y = Yh(:);

% INDICES
% -------------------------------------------------------------------------
Idx = reshape(1:((N + 2)^2), [N+2 N+2]);

Idx_xy           = Idx((2:N+1)    , (2:N+1)    );                           Idx_xy = Idx_xy(:); % (i  , j  )
Idx_xM(1:N, 1:N) = Idx((2:N+1) - 1, (2:N+1)    );                           Idx_xM = Idx_xM(:); % (i-1, j  )
Idx_xP(1:N, 1:N) = Idx((2:N+1) + 1, (2:N+1)    );                           Idx_xP = Idx_xP(:); % (i+1, j  )
Idx_yM(1:N, 1:N) = Idx((2:N+1)    , (2:N+1) - 1);                           Idx_yM = Idx_yM(:); % (i  , j-1)
Idx_yP(1:N, 1:N) = Idx((2:N+1)    , (2:N+1) + 1);                           Idx_yP = Idx_yP(:); % (i  , j+1)

Idx_Bd = Idx( Xh < h/2 | Yh < h/2 | Xh > L-h/2 | Yh > L-h/2 );              Idx_Bd = Idx_Bd(:); % index boundary points

%                                                   DIFFERENTIATION MATRIX
% -------------------------------------------------------------------------
% Sparse row indices
II = [Idx_xy;
      Idx_xy;
      Idx_xy;
      Idx_xy;
      Idx_xy;
      Idx_Bd];
  
% Sparse column indices
JJ = [Idx_xy;
      Idx_xM;
      Idx_xP;
      Idx_yM;
      Idx_yP;
      Idx_Bd]; % sparse column j

% Sparse values
AA = [ones(N^2, 1) * +4 / h^2;... % u(x,y)
      ones(N^2, 1) * -1 / h^2;... % u(x-1,y)
      ones(N^2, 1) * -1 / h^2;... % u(x+1,y)
      ones(N^2, 1) * -1 / h^2;... % u(x,y-1)
      ones(N^2, 1) * -1 / h^2;... % u(x,y+1)
      ones(length(Idx_Bd),1) ];
  
% Create sparse matrix
Lh = sparse(II, JJ, AA);