function [Xh, Lh, Fh] = a03ex03getBVP(P)
% ------------------------------------------------------------------------\
% Assignment 3, Exercise 2b                                               |
%                                                             submitted by|
%                                                                         |
%                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
%                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
%        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
%                                                                         |
%                                                        in  MATLAB R2014a|
% ------------------------------------------------------------------------/
%
%                                                                 Solution
% -------------------------------------------------------------------------
% Set initial boundary value paramters
% ------------------------------------
a     =  1;
b     = -4;
c     =  1;
alpha =  1;
beta  =  2;

% Set point resolution
n = 2.^P - 1;

% Set number of points
numberPoints = n + 1;

% Caclulate the step size
h = 1 / numberPoints;

% One dimensional coordinate vector
Xh = linspace(0, 1, numberPoints)' * h;

% Get F.D. stencils for m = 1
S = a03ex01getstencil(1);

% Copy the first three difference stencils
% ----------------------------------------
diff0thOrd = S(1,:);
diff1stOrd = S(2,:);
diff2ndOrd = S(3,:);

% Generate diagonal of the zeroth derivative matrix
L0 = c * spdiags(repmat(diff0thOrd, numberPoints, 1),... Row data of the diagonal matrix
                 -1:1                               ,... Position index of the diagonals
                 numberPoints                       ,... Number of matrix rows
                 numberPoints                       ); % Number of matrix columns

% Generate diagonal of the first derivative matrix
L1 =  b / h   * spdiags(repmat(diff1stOrd, numberPoints, 1),... Row data of the diagonal matrix
                        -1:1                               ,... Position index of the diagonals
                        numberPoints                       ,... Number of matrix rows
                        numberPoints                       ); % Number of matrix columns

% Generate diagonal of the second derivative matrix
L2 = -a / h^2 * spdiags(repmat(diff2ndOrd, numberPoints, 1),... Row data of the diagonal matrix
                        -1:1                               ,... Position index of the diagonals
                        numberPoints                       ,... Number of matrix rows
                        numberPoints                       ); % Number of matrix columns

% Generate the difference matrix
Lh = L1 + L2 + L0;

% Calculate the right hand side
Fh = Lh * ( 1 + 4 * Xh.^2 - 3 * Xh.^3 );

% Apply boundary conditions
% -------------------------
%{
Fh(1)   = Fh(1)   + alpha * (a / h^2 + b / (2*h));
Fh(end) = Fh(end) + beta  * (a / h^2 - b / (2*h));
%}