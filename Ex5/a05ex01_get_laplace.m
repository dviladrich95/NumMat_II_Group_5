function [Lh, XX, YY, BarOmega_h, Omega_h, Gamma_h] = a05ex01_get_laplace(L,N,IS_IN_DOMAIN)
% ------------------------------------------------------------------------\
% Assignment 5, Exercise 1                                                |
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

% generate coordinates and grid spacing
[XX,YY] = meshgrid(linspace(0,L,N));
XX = XX';
YY = YY';
h  = L / (N-1);

%                                                         INITIAL MATRICES
% -------------------------------------------------------------------------
BarOmega_h      = IS_IN_DOMAIN(XX,YY,L);
Omega_h         = logical(zeros(N,N));
Gamma_h         = logical(zeros(N,N));
Idx             = zeros(N,N);
Idx_xM          = zeros(N,N);
Idx_xP          = zeros(N,N);
Idx_yM          = zeros(N,N);
Idx_yP          = zeros(N,N);
Idx_Bd          = zeros(N,N);
Idx(BarOmega_h) = (1 : sum(BarOmega_h(:)));

%                                                      BOUNDARY ALLOCATION
% -------------------------------------------------------------------------
%{-
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
%}

%                           (is being kept as back-up) BOUNDARY ALLOCATION
% -------------------------------------------------------------------------
%{
%  Loop over all matrix columns
for i = 1 : size(BarOmega_h,1)
    % Loop over all matrix rows
    for j = 1 : size(BarOmega_h,2)
        % Condition for not in domain
        if ~BarOmega_h(i,j)
            Gamma_h(i, j) = false;
            Omega_h(i, j) = false;
            
        % Condition for gamma
        elseif (~BarOmega_h(i,j + 1) || ~BarOmega_h(i + 1,j) ||...
                ~BarOmega_h(i,j - 1) || ~BarOmega_h(i - 1,j))
            
            Gamma_h(i, j) = true;
            Omega_h(i, j) = false;
            Idx_Bd (i, j) = Idx(i, j);
            
        % Condition for omega
        else
            Gamma_h(i, j) = false;
            Omega_h(i, j) = true;
            
            Idx_xM(i, j) = Idx(i - 1, j    );
            Idx_xP(i, j) = Idx(i + 1, j    );
            Idx_yM(i, j) = Idx(i    , j - 1);
            Idx_yP(i, j) = Idx(i    , j + 1);
            
        end % of domain conditions
        
    end % of loop over all matrix rows
    
end % of loop over all matrix columns
%}

% Plot domain
% -----------
%{
figure(1)
hold on
plot(XX(~BarOmega_h),YY(~BarOmega_h),'o', 'MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
plot(XX(Omega_h)    ,YY(Omega_h)    ,'o', 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
plot(XX(Gamma_h)    ,YY(Gamma_h)    ,'o', 'MarkerFaceColor','b', 'MarkerEdgeColor', 'b');
hold off
figScaleFac = 0.75;
grid on
axis equal
xlabel('x', 'FontSize', 15)
ylabel('y', 'FontSize', 15)
set(gca, 'FontSize', 15)
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
%}

%                                                   DIFFERENTIATION MATRIX
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