function [Ah, Mh, Fh] = a07ex01getPDE(XH, F, CONST, THETA, TAU, FLAG)
% ------------------------------------------------------------------------\
% Assignment 7, Exercise 1c                                               |
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
% Computation of the right hand side vector
Fh = F(XH); Fh = Fh';

% Computation of the Laplace operator
Lh = getLaplace(length(XH),... Number of points
               XH(2)-XH(1),... Spacing
                  CONST(1),... a
                  CONST(2),... b
                         0,... c
                      FLAG); % Difference stencil
                     
% Computation of the Ah matrix               
Ah = eye(length(Lh)) + TAU * (THETA - 0) * Lh;

% Computation of the Mh matrix
Mh = eye(length(Lh)) + TAU * (THETA - 1) * Lh;

%                                                           LAPLACE MATRIX
% -------------------------------------------------------------------------
function Lh = getLaplace(N, H, A, B, C, FLAG)
% Difference conditions for u'
switch FLAG
    % Condition for forward difference
    case '+'
        % Difference stencils of the first two derivatives
        S = [ 1,  0, 0;
             -1,  1, 0;
              1, -2, 1];

        % Diagonal positions of the stencil coefficients
        diagCol = [0, 1, 2];
        
    % Condition for backward difference
    case '-'
        % Difference stencils of the first two derivatives
        S = [ 1,  0, 0;
              1, -1, 0;
              1, -2, 1];
  
        % Diagonal positions of the stencil coefficients
        diagCol = [0, -1, -2];
        
    % Condition for central difference
    case '0'
        % Difference stencils of the first two derivatives
        S = [   0,  1,    0; 
             -1/2,  0,  1/2;
                1, -2,    1];

        % Diagonal positions of the stencil coefficients
        diagCol = [-1, 0, 1];
        
end % of difference conditions

% Difference stencil for u''
S2 = [1, -2,    1];

% Diagonal positions of the stencil coefficients for u''
diagCol2 = [-1 0 1];

% Generate diagonal of the zeroth derivative matrix
L0 = C * spdiags(repmat(S(1,:), N, 1),... Row data of the diagonal matrix
                 diagCol             ,... Position index of the diagonals
                 N                   ,... Number of matrix rows
                 N                   ); % Number of matrix columns

% Generate diagonal of the first derivative matrix
L1 =  B / H   * spdiags(repmat(S(2,:), N, 1),... Row data of the diagonal matrix
                        diagCol             ,... Position index of the diagonals
                        N                   ,... Number of matrix rows
                        N                   ); % Number of matrix columns

% Generate diagonal of the second derivative matrix
L2 = -A / H^2 * spdiags(repmat(S2, N, 1),... Row data of the diagonal matrix
                        diagCol2        ,... Position index of the diagonals
                        N               ,... Number of matrix rows
                        N               ); % Number of matrix columns

% Generate the difference matrix
Lh = L0 + L1 + L2;

% Adjust the diagonal elements on boundary points
% -----------------------------------------------
Lh(  1, :) = 0;
Lh(end, :) = 0;