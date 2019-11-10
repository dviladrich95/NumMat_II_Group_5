function S = a03ex01getstencil(M)
% ------------------------------------------------------------------------\
% Assignment 3, Exercise 1b                                               |
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
% Generate matrix for point index
J = meshgrid([- M : 1 : M])';

% Generate matrix for derivative order
K = meshgrid([0 : 1 : 2*M]);

% Generate the stencil matrix
S = inv(J.^K ./  factorial(K));