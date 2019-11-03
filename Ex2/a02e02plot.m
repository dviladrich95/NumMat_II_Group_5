function a02e02plot
% ------------------------------------------------------------------------\
% Assignment 2, Exercise 2                                                |
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
% Set number of points
numberPoints = 50;

% Set time vector
TIME = [0, 1E-3, 1E-2, 1E-1, 1];

% Set X vector
XPOINT = linspace(0, 1, numberPoints)';

% Cases:      a | c | mode
InputCase = [ 1,  0,  1;
              1,  1,  1;
              1,  0,  2;
              1,  1,  2];

% Loop over all cases
for caseIdx = 1 : size(InputCase, 1)
    % Open and configure new figure
    % -----------------------------
    figure(caseIdx)
    hold all
    grid on
    
    % Loop over all time values
    for timeIdx = 1 : length(TIME)
        % Solve the heat equation for current parameters
        U = a02e02solveheat(InputCase(caseIdx, 3),... mode
                            XPOINT               ,... Vector of X points
                            TIME(timeIdx)        ,... current time
                            InputCase(caseIdx, 1),... a
                            InputCase(caseIdx, 2)); % c
        
        % Plot the temperature over time
        plot(XPOINT, U, 'LineWidth', 2);
                                    
    end % of loop over all time values
    
    % Plot settings
    % -------------
    title(['Mode:' num2str(InputCase(caseIdx, 3)) ', ' ...
           'a='    num2str(InputCase(caseIdx, 1)) ', ' ...
           'c='   num2str(InputCase(caseIdx, 2))      ], 'FontSize', 15)
    legend('t = 0', 't = 0.001', 't = 0.01', 't = 0.1', 't = 1')
    xlabel('x', 'FontSize', 15)
    ylabel('u(x,t)', 'FontSize', 15)
    set(gcf, 'Color', 'white')
    set(gca, 'FontSize', 15)
    
end % of loop over all cases

%{
InputCase  = [ -1,  0,  1;
               -1,  0,  2];
%}