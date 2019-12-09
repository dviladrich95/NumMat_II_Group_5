function a06ex01_solveEVPc
% ------------------------------------------------------------------------\
% Assignment 6, Exercise 1c                                               |
%                                                             submitted by|
%                                                                         |
%                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
%                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
%                   Riccardo Parise | 412524 | Scientific Computing, M.Sc.|
%        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
%                                                                         |
%                                                        in  MATLAB R2014a|
% ------------------------------------------------------------------------/
close all
clear
clc
%                                                                 Solution
% -------------------------------------------------------------------------
% Control parameters
% ------------------
N =  2.^([4, 6, 9]) - 1;
l =  2.2;
r =    1;
nEV     = 6;
plotIdx = 0;

% Mesh grid
% ---------
[XX,YY] = meshgrid(linspace(0,l,N(end)));
XX = XX' - l/2;
YY = YY' - l/2;
RR = sqrt(XX.^2 + YY.^2);
PHIPHI = atan2(YY,XX);
VV = nan(N(end),N(end));

% Get analytical eigen values
[EV(:,1), Nu] = a06ex01_solveEVPb(1);

% Loop over all meshes
for n = N
    % Get reduced laplace matrix for Dirichlet boundary conditions
    [Lh, ActiveDOFs, BarOmega_h] = getLaplaceDirichlet(l, n, r);
    
    % Calculate numerical eigen values
    EV(:, end + 1) = eigs(Lh(ActiveDOFs, ActiveDOFs), nEV, 'sm');
    
end % of loop over all meshes


%                                                           Print solutions
% -------------------------------------------------------------------------
% Print header of the table
% -------------------------
fprintf('/---------------------------------------------------------------------\\\n')
fprintf('|\t       N = %d |\t      N = %d |\t     N = %d |\t  Analytical |\t# |\n', N(1),...
                                                                               N(2),...
                                                                               N(3))
fprintf('| --------------------------------------------------------------------|\n')

% Loop over all eigen values
for evCount = 1 : nEV
    % Print eigenvalues
    fprintf('|\t%d |\t%d |\t%d |\t%d |\t%i |\n', EV(evCount,2),...
                                                 EV(evCount,3),...
                                                 EV(evCount,4),...
                                                 EV(evCount,1),...
                                                 nEV - evCount +1);
    
end % of loop over all eigenvalues

% End table
fprintf('\\---------------------------------------------------------------------/\n')

% Eigenfunctions
for evIdx = [1, 2, 4, 6];
    plotIdx = plotIdx + 1;
    subplot(2,2, plotIdx)
    VV(BarOmega_h) = cos(Nu(evIdx) * PHIPHI(BarOmega_h)) .* besselj(Nu(evIdx), RR(BarOmega_h) * sqrt(EV(evIdx,1)));
    surf(XX,YY,VV, 'EdgeAlpha', 0)
    title(['$\lambda_' num2str(evIdx) ' = ' num2str(EV(7 - evIdx,1)) '$'], 'FontSize', 13, 'Interpreter','latex')
    axis tight equal
    view([0,90])
    xlabel('x', 'FontSize', 15)
    ylabel('y', 'FontSize', 15)
    colorbar
    set(gca, 'FontSize', 15)
end

figScaleFac = 0.75;
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
%export_fig('../Documentation/Figures/a06ex01c.png')

%                                         FUNCTION FOR LAPLACE - DIRICHLET
% -------------------------------------------------------------------------
function [Lh, ActiveDOFs, BarOmega_h] = getLaplaceDirichlet(L, N, R)
% Generate coordinates and grid spacing
% -------------------------------------
[XX,YY] = meshgrid(linspace(0,L,N));
XX = XX';
YY = YY';
h  = L / (N-1);

%                                                         Initial Matrices
% -------------------------------------------------------------------------
BarOmega_h      = ((XX - L/2).^2 + (YY - L/2).^2) <= R;
Omega_h         = logical(zeros(N,N));
Gamma_h         = logical(zeros(N,N));
Idx             = zeros(N,N);
Idx_xM          = zeros(N,N);
Idx_xP          = zeros(N,N);
Idx_yM          = zeros(N,N);
Idx_yP          = zeros(N,N);
Idx_Bd          = zeros(N,N);
Idx(BarOmega_h) = (1 : sum(BarOmega_h(:)));

%                                                      Boundary Allocation
% -------------------------------------------------------------------------
% Set the slicing range
SliceRange = 2 : N-1;

% Find points on gamma_h
Gamma_h(SliceRange,SliceRange) = BarOmega_h(SliceRange    ,     SliceRange) &...
                                 (BarOmega_h(SliceRange + 1, SliceRange    ) == 0 |...
                                  BarOmega_h(SliceRange - 1, SliceRange    ) == 0 |...
                                  BarOmega_h(SliceRange    , SliceRange + 1) == 0 |...
                                  BarOmega_h(SliceRange    , SliceRange - 1) == 0);

% Find points in omega_h
Omega_h(SliceRange,SliceRange) = BarOmega_h(SliceRange, SliceRange) &...
                                 ~Gamma_h(SliceRange, SliceRange);
                         
% Get neigbour indices
% --------------------
Idx_xM(SliceRange, SliceRange) = Idx(SliceRange - 1, SliceRange    );
Idx_xP(SliceRange, SliceRange) = Idx(SliceRange + 1, SliceRange    );
Idx_yM(SliceRange, SliceRange) = Idx(SliceRange    , SliceRange - 1);
Idx_yP(SliceRange, SliceRange) = Idx(SliceRange    , SliceRange + 1);

%                                                   Differentiation Matrix
% -------------------------------------------------------------------------
% Sparse row indices
II = [Idx(Omega_h);
      Idx(Omega_h);
      Idx(Omega_h);
      Idx(Omega_h);
      Idx(Omega_h);
      Idx(Gamma_h)];
  
% Sparse column indices
JJ = [Idx(Omega_h)   ;
      Idx_xM(Omega_h);
      Idx_xP(Omega_h);
      Idx_yM(Omega_h);
      Idx_yP(Omega_h);
      Idx(Gamma_h)   ];

% Sizes of the domain and boundary vectors
% ----------------------------------------
k1 = numel(find(Omega_h));
k2 = numel(find(Gamma_h));

% Sparse values
AA = [ones(k1, 1) * +4 / h^2;... % u(x,y)
      ones(k1, 1) * -1 / h^2;... % u(x-1,y)
      ones(k1, 1) * -1 / h^2;... % u(x+1,y)
      ones(k1, 1) * -1 / h^2;... % u(x,y-1)
      ones(k1, 1) * -1 / h^2;... % u(x,y+1)
      ones(k2, 1) ];

% Create sparse matrix
Lh = sparse(II, JJ, AA);

% Get indices for active DOFs
ActiveDOFs = Idx(Omega_h);