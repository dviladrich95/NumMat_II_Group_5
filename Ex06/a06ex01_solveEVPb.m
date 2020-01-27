function [EigVal, Nu] = a06ex01_solveEVPb(R)
% ------------------------------------------------------------------------\
% Assignment 6, Exercise 1b                                               |
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
% Matrix of the k-th Bessel zeros
Bessel = [2.40483, 5.52008,  8.65373, 11.79153, 14.93092; % n = 0
          3.83171, 7.01559, 10.17347, 13.32369, 16.47063; % n = 1
          3.83171, 7.01559, 10.17347, 13.32369, 16.47063; % n = 1 - duplicated
          5.13562, 8.41724, 11.61984, 14.79595, 17.95982; % n = 2
          5.13562, 8.41724, 11.61984, 14.79595, 17.95982; % n = 2 - duplicated
          6.38016, 9.76102, 13.01520, 16.22347, 19.40942; % n = 3
          6.38016, 9.76102, 13.01520, 16.22347, 19.40942];% n = 3 - duplicated
   
% Discretisized eigen values of the given disc problem for given radius
EigVal = sort(Bessel(:)).^2 / R^2;

% Crop first six eigenvalues
EigVal = flipud(EigVal(1:6));

% Nu values of the first six eigen values
Nu = fliplr([0, 1, 1, 2, 2, 0]);