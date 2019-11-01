function a02e01plot
% ------------------------------------------------------------------------\
% Assignment 2, Exercise 1                                                |
%                                                             submitted by|
%                                                                         |
%                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
%                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
%        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
%                                                                         |
%                                                        in  MATLAB R2014a|
% ------------------------------------------------------------------------/
close all
clear
clc
%                                                                 Solution
% -------------------------------------------------------------------------
% Set the dimension of data points
numberPoints = 150; 

% Generate v(r, phi) data
% -----------------------
R   = linspace(1,    2, numberPoints);
PHI = linspace(0, 2*pi, numberPoints);

% Generate the mesh grid of 
[RR, PHIPHI] = meshgrid(R, PHI);

% Conversion of polar coordinates to cartesian coordinates
% --------------------------------------------------------
XX = RR .* cos(PHIPHI);
YY = RR .* sin(PHIPHI);

% Surface mesh of the disc with solution
surfID = surf(XX                    ,...
              YY                    ,...
              a02e01solution(XX, YY));

% Surface plot settings
% ---------------------
axis equal
xlabel('X')
ylabel('Y')
zlabel('U(x,y)')
set(surfID, 'EdgeAlpha', .2)