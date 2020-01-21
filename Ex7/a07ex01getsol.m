function [U] = a07ex01getsol(X_H, T)
% ------------------------------------------------------------------------\
% Assignment 7, Exercise 1d                                               |
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
% Fourrier row indices
K = 1 : 1 : 300;

% Loop over all Fourrier rows
for k = K
    % Calculate inhomogenous problem
    U0(k, :) = (-1)^(k-1) / (k * pi) * sin(k * 2 * pi * X_H) * exp(-(k*2*pi)^2 * T);
    
end % of loop over Fourrier rows

% Calculate U by summing all Fourrier rows and subtract from the homogenous part
U = X_H - sum(U0, 1); U = U';