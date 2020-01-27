function [Lh, Xh, Yh, Idx_Bd] = a04ex01_Lh5(L, N)
% ------------------------------------------------------------------------\
% Assignment 4, Exercise 1c                                               |
%                                                             submitted by|
%                                                                         |
%                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
%                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
%                   Riccardo Parise | 412524 | Scientific Computing, M.Sc.|
%        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
%                                                                         |
%                                                        in  MATLAB R2014a|
% ------------------------------------------------------------------------/
%
%                                                                 Solution
% -------------------------------------------------------------------------
% Caclulate the step size
h = L / (N + 1);

% Generate meshgrid of the domain
[Xh, Yh] = meshgrid(linspace(0, L, N + 2));

% Transpose the mesh grid
Xh = Xh'; Yh = Yh';

%                                                                 INDEXING
% -------------------------------------------------------------------------
% Generate index data of the domain
Idx = reshape(1:((N + 2)^2), [N+2 N+2]);

% Generate index data for points that are not on the boundary
Idx_xy           = Idx((2:N+1)    , (2:N+1)    );                           Idx_xy = Idx_xy(:); % (i  , j  )

% % Generate index data for five point stencil negihbors
% ------------------------------------------------------
Idx_xM(1:N, 1:N) = Idx((2:N+1) - 1, (2:N+1)    );                           Idx_xM = Idx_xM(:); % (i-1, j  )
Idx_xP(1:N, 1:N) = Idx((2:N+1) + 1, (2:N+1)    );                           Idx_xP = Idx_xP(:); % (i+1, j  )
Idx_yM(1:N, 1:N) = Idx((2:N+1)    , (2:N+1) - 1);                           Idx_yM = Idx_yM(:); % (i  , j-1)
Idx_yP(1:N, 1:N) = Idx((2:N+1)    , (2:N+1) + 1);                           Idx_yP = Idx_yP(:); % (i  , j+1)

% Generate index data boundary points
Idx_Bd = Idx( Xh < h/2 | Yh < h/2 | Xh > L-h/2 | Yh > L-h/2 );              Idx_Bd = Idx_Bd(:); % Boundary indices

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