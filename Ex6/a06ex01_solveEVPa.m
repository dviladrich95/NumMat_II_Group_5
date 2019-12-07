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
m = 2; % Number of the neighbour points in one direction. WARNING: Taking m = 2 (5 point stencil) leads to flip of the numerical eigenlavues of the printed lists.
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
[EVec1, EV1_Num] = eigs(Lh(ActiveDOF, ActiveDOF), nEV, 'sm');

% ANalytical solution
EV1_Exact = ([1:n]*pi).^2 + 1;

%                                                  a.2) Homogenous Neumann
% -------------------------------------------------------------------------
% Get laplace for Neumann boundary conditions
Lh = getLaplaceNeumann(Lh, ActiveDOF, h);
                   
% Calculate the first 20 eigenvalues
[EVec2, EV2_Num] = eigs(Lh, nEV, 'sm');

% Analytical solution
EV2_Exact = ([0:n]*pi).^2 + 1;

%                                                            a.3) Periodic
% -------------------------------------------------------------------------""
% Get laplace for periodic boundary conditions
Lh = getLaplacePeriodic(n, m, h, a, b, c);

% Calculate the first 20 eigenvalues
[EVec3, EV3_Num] = eigs(Lh, 2*nEV, 'sm');

% Anaytical solution
EV3_Exact = (2 * [0:n-1] * pi).^2 + 1;

%                                                        Plot eigenvectors
% -------------------------------------------------------------------------
%{-
% Some text paramters
% -------------------
labelFont = 20;
axisFont  = 15;

figure(1)
subplot(1,3,1)
plot(h*[1:n],   EVec1(:,end-1)', '-r',...
     h*[1:n], sin(2*pi*[1:n]*h),'--k',...
     'LineWidth', 2)
 
grid on
title('Dirichlet', 'FontSize', labelFont)
xlabel('x', 'FontSize', labelFont)
ylabel('u', 'FontSize', labelFont)

subplot(1,3,2)
plot(h*[1:n],   EVec2(:,end-2)', '-r',...
     h*[1:n], cos(2*pi*[1:n]*h),'--k',...
     'LineWidth', 2)

grid on
title('Neumann', 'FontSize', labelFont)
xlabel('x', 'FontSize', labelFont)
ylabel('u', 'FontSize', labelFont)
 
subplot(1,3,3)
plot(h*[0:n],   EVec3(:,end-4)', '-r',...
     h*[1:n], cos(2*2*pi*[1:n]*h),'--k',...
     'LineWidth', 2)
 
grid on
title('Periodic', 'FontSize', labelFont)
xlabel('x', 'FontSize', labelFont)
ylabel('u', 'FontSize', labelFont)
 
figScaleFac = 0.75;
set(gca, 'FontSize', 15)
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
%export_fig('../Documentation/Figures/a06ex01EV2.png')
%}

%                                                        Print eigenvalues
% -------------------------------------------------------------------------
% Print header of the table
% -------------------------
fprintf('/-------------------------------------------------------------------------------------------------------\\\n')
fprintf('|\t # |\t  Dir. Num. |\t Dir. Exact |\t  Neu. Num. |\t Neu. Exact |\tPrdic. Num. |\tPrdic. Exact|\n')
fprintf('|-------------------------------------------------------------------------------------------------------|\n')

% Loop over all eigen values
for evCount = 1 : nEV
    % Print eigenvalues
    fprintf('|\t%2i |\t%e|\t%e|\t%e|\t%e|\t%e|\t%e|\n', nEV - evCount + 1    ,...
                                              EV1_Num(evCount, evCount),...
                                              EV1_Exact(nEV - evCount + 1),...
                                              EV2_Num(evCount, evCount),...
                                              EV2_Exact(nEV - evCount + 1),...
                                              EV3_Num(2*evCount, 2*evCount),...
                                              EV3_Exact(nEV - evCount + 1));
    
end % of loop over all eigenvalues

% End table
fprintf('\\-------------------------------------------------------------------------------------------------------/\n')

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
function [Lh] = getLaplaceNeumann(Lh, ACTIVE_DOF, H)
% Adjust Lh for Neumann
% ---------------------
Lh = Lh(ACTIVE_DOF, ACTIVE_DOF);

Lh(1,1) = 1/H^2; Lh(1,2) = -1/H^2;
Lh(end,end-1) = -1/H^2; Lh(end,end) = 1/H^2;


%                                           FUNCTION FOR LAPLCE - PERIODIC
% -------------------------------------------------------------------------
function [Lh] = getLaplacePeriodic(N, M, H, A, B, C)
% Set number of points
numberPoints = N+1;

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