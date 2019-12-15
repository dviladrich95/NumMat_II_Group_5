% ------------------------------------------------------------------------\
% Assignment 7, Exercise 1d)                                              |
%                                                             submitted by|
%                                                                         |
%                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
%                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
%                   Riccardo Parise | 412524 | Scientific Computing, M.Sc.|
%        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
%                                                                         |
%                                                        in  MATLAB R2014a|
% ------------------------------------------------------------------------/
close all
clear
clc
%                                                                 Solution
% -------------------------------------------------------------------------
% PDE Parameters
% --------------
f = @(X) zeros(1, length(X));
Theta = [0, 0.5, 1];
Const = [1, 0];
N     = [10, 100, 1000];
M     = [50, 500, 5000];
t     = 0.01;

% Set initial figure index
figIdx = 0;

% Loop over all theta schemes
for theta = Theta
    % Get to next figure
    figIdx = figIdx + 1;
    
    % Open the figure
    figure(figIdx)

    % Plot windows settings
    % ---------------------
    figScaleFac = 1;
    set(gcf, 'Color'            , 'white'                                          ,...
             'PaperSize'        , [34, 34]                                         ,...
             'PaperPositionMode', 'auto'                                           ,...
             'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac * 1])
    
    % Reset n counter
    nIdx = 0;
    
    % Loop over all grid resolutions
    for n = N
        % Expand n for boundary points
        n = n + 2;
                
        % Get to n
        nIdx = nIdx + 1;
                
        % Reset m counter
        mIdx = 0;
        
        % Loop over all time step resolutions
        for m = M
            % Get to next m
            mIdx = mIdx + 1;
            
            % Get to next subplot
            subplot(3,3, (nIdx-1) * 3 + mIdx)
                        
            % Generate 1D arrays
            % ------------------
            Xh = linspace(0, 1, n);
            Th = linspace(1E-5, t, m);

            % Get PDE
            [Ah, Mh, Fh] = a07ex01getPDE(Xh, f, Const, theta, Th(2)-Th(1), '0');

            % Solve PDE as numerical and analtical
            [Um, err] = a07ex01d(n, m, t, Ah, Mh, Fh);
                                          
            % Stability measure ( must be greater than 1, II.64)
            stab = 1 - 2 * (1 - theta) * (Th(2)-Th(1)) / (Xh(2)-Xh(1))^2;
            
            % Calculate and print maximum norm
            fprintf('Theta: %d, n: %i, m: %i, 1-2r = %d, Error: %d\n', theta, n, m, stab, err);
        
        end % of loop over all time step resolutions

    end % of loop over all grid resolutions
    
    % Remark
    disp('THE REULTS WERE DISCUSSED IN THE ASSINGMENT PDF!')

    
    % Save figure
    %export_fig(['../Documentation/Figures/a07ex01d_Theta' num2str(figIdx) '.png'])
    
end % of loop over all theta schemes