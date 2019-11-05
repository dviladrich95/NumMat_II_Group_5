function [Xh, Lh] = a03ex01getLaplace(M,N)
% ------------------------------------------------------------------------\
% Assignment 3, Exercise 1c                                               |
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
% Set number of points
numberPoints = N + 1;

% Step size
h = 1 / numberPoints;

% One dimensional coordinate vector
Xh = linspace(0, N*h, numberPoints);

% Stencil matrix
S = a03ex01getstencil(M);

% Stencil constants for the second order derivative
diffOrd = S(3,:);

% Generate the main diagonal matrix
LDiag = spdiags(repmat(diffOrd, numberPoints, 1),... Row data of the diagonal matrix
                -M:M                            ,... Position index of the diagonals
                numberPoints                    ,... Number of matrix rows
                numberPoints                    ); % Number of matrix columns
         
% Upper triangle matrix
LUpper = spdiags(repmat(diffOrd,numberPoints, 1),... Row data of the diagonal matrix
                 (numberPoints-M):numberPoints  ,... Position index of the diagonals
                 numberPoints                   ,... Number of matrix rows
                 numberPoints                   ); % Number of matrix columns

% Lower triangle matrix
LLower = spdiags(repmat(fliplr(diffOrd),numberPoints, 1)       ,... Row data of the diagonal matrix
                 fliplr((-numberPoints + 1):(M - numberPoints)),... Position index of the diagonals
                 numberPoints                                  ,... Number of matrix rows
                 numberPoints                                  ); % Number of matrix columns
         
         
% Summation of sub matrices for second order derivative
Lh = (LDiag + LUpper + LLower) / h^2;