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

% Vectors for difference resolution
% ---------------------------------
N = 2.^(2:15) - 1;
H = 1 ./ (N + 1);

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
    plot(H, DiffErr, 'LineWidth', 2)
    
    % Add legend data for current m
    legendData{end +1} = ['m = ' num2str(m)];
        
end % of loop over all M

% Plot settings
% -------------
figScaleFac = 0.75;
grid on
legend(legendData)
xlabel('h', 'FontSize', 15)
ylabel('Error', 'FontSize', 15)
set(gca, 'FontSize', 15       ,...
         'XScale'  , 'log'    ,...
         'YScale'  , 'log'    ,...
         'XDir'    , 'reverse')
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])

%{
    addpath('C:\exportFig\')
    export_fig('../Documentation/Figures/a03e01DiffError')
%}