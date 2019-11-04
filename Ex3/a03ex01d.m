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
close all
clear
clc
%                                                                 Solution
% -------------------------------------------------------------------------
% Set initial figure data
% -----------------------
figure(1)
hold all
legendData = {};

% Number of neighbour points
M = 1:4;

% Vector for difference resolution
N = 2.^(2:15) - 1;

% Loop over all M
for m = M
    % Set initial error vector
    DiffErr = [];
    
    % Loop over all P
    for n = N        
        % Get X and laplace matrix for the second order derivative
        [Xh, Lh] = a03ex01getLaplace(m, n);
                
        % Calculate second order derivative analyticaly
        U2d = - (2*pi)^2 * sin(2*pi*Xh)';
        
        % Calculate the infinity norm of the error between numerical and  analytical derivatives
        DiffErr(end + 1) = max(abs(Lh * sin(2*pi*Xh)' - U2d));
        
    end % of loop over all p
    
    % Plot error over difference resolution
    plot(N,DiffErr)
    
    % Add legend data for current m
    legendData{end +1} = ['M = ' num2str(m)];
        
end % of loop over all M

% Plot settings
% -------------
legend(legendData)