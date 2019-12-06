function a06ex01_solveEVPa
% ------------------------------------------------------------------------\
% Assignment 6, Exercise 1a                                               |
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
n = 200;
m = 1;
a = 1;
b = 0;
c = 1;
nEV = 20;

% Active indices
ActiveDOF = 2 : n+1;

% Unit matrix for eigenvalue calculations
M = speye(n+2);

% Step size
h = 1 / (n+1);

%                                                a.1) Homogenous Dirichlet
% -------------------------------------------------------------------------
% Get laplace for Dirichlet boundary conditions
Lh = getLaplaceDirichlet(n, m, h, a, b, c);
                  
% Calculate the first 20 numerical eigenvalues
[EVec1, EV1] = eigs(Lh(ActiveDOF, ActiveDOF), M(ActiveDOF, ActiveDOF), nEV, 'sm');

%                                                  a.2) Homogenous Neumann
% -------------------------------------------------------------------------
% Get laplace for Neumann boundary conditions
Lh = getLaplaceNeumann(Lh);
                   
% Calculate the first 20 eigenvalues
[EVec2, EV2] = eigs(Lh(ActiveDOF, ActiveDOF), M(ActiveDOF, ActiveDOF), nEV, 'sm');

%                                                            a.3) Periodic
% -------------------------------------------------------------------------
% Get laplace for periodic boundary conditions
Lh = getLaplacePeriodic(n, m, h, a, b, c);

% Calculate the first 20 eigenvalues
[EVec3, EV3] = eigs(Lh(ActiveDOF, ActiveDOF), M(ActiveDOF, ActiveDOF), nEV, 'sm');

%                                                      Analytical Solution
% -------------------------------------------------------------------------
% Calculate the first 20 analytical eigenvalues
EV_ana= 4 / (h^2) * sin([fliplr(1:nEV)]' * h *  pi / 2).^2 + 1;

%                                                        Plot eigenvectors
% -------------------------------------------------------------------------
% 2nd Eigenvalue
% --------------
figure(1)
plot(h*ActiveDOF,        EVec1(:,end-1), '-r',...
     h*ActiveDOF,        EVec2(:,end-1), '-b',...
     h*ActiveDOF,        EVec3(:,end-1), '-g',...
     h*ActiveDOF, sin(2*pi*h*ActiveDOF),'--k',...
     'LineWidth', 2)
grid on
figScaleFac = 0.75;
xlabel('x', 'FontSize', 25)
ylabel('u', 'FontSize', 25)
set(gca, 'FontSize', 15)
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
%export_fig('../Documentation/Figures/a06ex01EV2.png')
 
% 20th Eigenvalue
% ---------------
figure(2)
plot(h*ActiveDOF,             EVec1(:,1), '-r',...
     h*ActiveDOF,             EVec2(:,1), '-b',...
     h*ActiveDOF,             EVec3(:,1), '-g',...
     h*ActiveDOF, sin(20*pi*h*ActiveDOF),'--k',...
     'LineWidth', 2)
grid on
figScaleFac = 0.75;
xlabel('x', 'FontSize', 25)
ylabel('u', 'FontSize', 25)
set(gca, 'FontSize', 15)
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
%export_fig('../Documentation/Figures/a06ex01EV20.png')

%                                                        Print eigenvalues
% -------------------------------------------------------------------------
% Print header of the table
% -------------------------
fprintf('/-----------------------------------------------------------------------\\\n')
fprintf('|\t # |\t Analytical |\t  Dirichlet |\t    Neumann |\t   Periodic |\n')
fprintf('|-----------------------------------------------------------------------|\n')

% Loop over all eigen values
for evCount = 1 : nEV
    % Print eigenvalues
    fprintf('|\t%2i |\t%d|\t%d|\t%d|\t%d|\n', nEV - evCount + 1    ,...
                                              EV_ana(evCount)      ,...
                                              EV1(evCount, evCount),...
                                              EV2(evCount, evCount),...
                                              EV3(evCount, evCount));
    
end % of loop over all eigenvalues

% End table
fprintf('\\-----------------------------------------------------------------------/\n')

%                                                FUNCTION FOR DIFF STENCIL
% -------------------------------------------------------------------------
function S = getStencil(M)
% Generate matrix for point index
J = meshgrid([- M : 1 : M])';

% Generate matrix for derivative order
K = meshgrid([0 : 1 : 2*M]);

% Generate the stencil matrix
S = inv(J.^K ./  factorial(K));

%                                         FUNCTION FOR LAPLACE - DIRICHLET
% -------------------------------------------------------------------------
function Lh = getLaplaceDirichlet(N, M, H, A, B, C)
% Set number of points
numberPoints = N + 2;

% Stencil matrix
S = getStencil(M);

% Generate diagonal of the zeroth derivative matrix
L0 = C * spdiags(repmat(S(1,:), numberPoints, 1),... Row data of the diagonal matrix
                 -M:M                           ,... Position index of the diagonals
                 numberPoints                   ,... Number of matrix rows
                 numberPoints                   ); % Number of matrix columns

% Generate diagonal of the first derivative matrix
L1 =  B / H   * spdiags(repmat(S(2,:), numberPoints, 1),... Row data of the diagonal matrix
                        -M:M                           ,... Position index of the diagonals
                        numberPoints                   ,... Number of matrix rows
                        numberPoints                   ); % Number of matrix columns

% Generate diagonal of the second derivative matrix
L2 = -A / H^2 * spdiags(repmat(S(3,:), numberPoints, 1),... Row data of the diagonal matrix
                        -M:M                           ,... Position index of the diagonals
                        numberPoints                   ,... Number of matrix rows
                        numberPoints                   ); % Number of matrix columns

% Generate the difference matrix
Lh = L0 + L1 + L2;

%                                           FUNCTION FOR LAPLACE - NEUMANN
% -------------------------------------------------------------------------
function [Lh] = getLaplaceNeumann(Lh)
% Adjust Lh for Neumann
% ---------------------
Lh          = [Lh; ones(1, length(Lh))];
Lh          = [Lh, ones(length(Lh), 1)];
Lh(end,end) = 0;

%                                           FUNCTION FOR LAPLCE - PERIODIC
% -------------------------------------------------------------------------
function [Lh] = getLaplacePeriodic(N, M, H, A, B, C)
% Set number of points
numberPoints = N + 1;

% Stencil matrix
S = getStencil(M);

% Generate zeroth derivative matrix
L0 = C * (spdiags(repmat(S(1,:), numberPoints, 1),... Row data of the diagonal matrix
                -M:M                            ,... Position index of the diagonals
                numberPoints                    ,... Number of matrix rows
                numberPoints                    )... Number of matrix columns
         ...     
        + spdiags(repmat(S(1,:),numberPoints, 1),... Row data of the diagonal matrix
         (numberPoints-M):numberPoints          ,... Position index of the diagonals
          numberPoints                          ,... Number of matrix rows
          numberPoints                          )... Number of matrix columns
         ...
        + spdiags(repmat(fliplr(S(1,:)),numberPoints, 1),... Row data of the diagonal matrix
         fliplr((-numberPoints + 1):(M - numberPoints)) ,... Position index of the diagonals
         numberPoints                                   ,... Number of matrix rows
         numberPoints                                   ));% Number of matrix columns

% Generate first derivative matrix
L1 = B / H * (spdiags(repmat(S(2,:), numberPoints, 1),... Row data of the diagonal matrix
                     -M:M                            ,... Position index of the diagonals
                     numberPoints                    ,... Number of matrix rows
                     numberPoints                    )... Number of matrix columns
             ...     
            + spdiags(repmat(S(2,:),numberPoints, 1),... Row data of the diagonal matrix
             (numberPoints-M):numberPoints          ,... Position index of the diagonals
              numberPoints                          ,... Number of matrix rows
              numberPoints                          )... Number of matrix columns
             ...
            + spdiags(repmat(fliplr(S(2,:)),numberPoints, 1),... Row data of the diagonal matrix
             fliplr((-numberPoints + 1):(M - numberPoints)) ,... Position index of the diagonals
             numberPoints                                   ,... Number of matrix rows
             numberPoints                                   ));% Number of matrix columns

% Generate first derivative matrix
L2 = -A / H^2 * (spdiags(repmat(S(3,:), numberPoints, 1),... Row data of the diagonal matrix
                      -M:M                            ,... Position index of the diagonals
                      numberPoints                    ,... Number of matrix rows
                      numberPoints                    )... Number of matrix columns
              ...     
             + spdiags(repmat(S(3,:),numberPoints, 1),... Row data of the diagonal matrix
              (numberPoints-M):numberPoints          ,... Position index of the diagonals
               numberPoints                          ,... Number of matrix rows
               numberPoints                          )... Number of matrix columns
              ...
             + spdiags(repmat(fliplr(S(3,:)),numberPoints, 1),... Row data of the diagonal matrix
              fliplr((-numberPoints + 1):(M - numberPoints)) ,... Position index of the diagonals
              numberPoints                                   ,... Number of matrix rows
              numberPoints                                   ));% Number of matrix columns
         
% Summation of sub matrices for second order derivative
Lh = L0 + L1 + L2;