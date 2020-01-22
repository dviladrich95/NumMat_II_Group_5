% ------------------------------------------------------------------------\
% Assignment 7, Exercise 1f)                                              |
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
f = @(X) [0, ones(1, length(X) - 2) ,0];
Theta = [ 0,  0.5,   1];
Flag  = {'0', '+', '-'};
Const = [1/100 1];
n     = 50 + 2;
m     = 100;
t     = 100;
t0    = 1E-5;

% Plot windows settings
% ---------------------
figScaleFac = 1;
set(gcf, 'Color'            , 'white'                                          ,...
         'PaperSize'        , [34, 34]                                         ,...
         'PaperPositionMode', 'auto'                                           ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac * 1])

% Reset theta index
thetaIdx = 0;
     
% Loop over all theta schemes
for theta = Theta
    thetaIdx = thetaIdx + 1;
    
    flagIdx = 0;
    % Loop over all difference stencils
    for flag = Flag
        flagIdx = flagIdx + 1;
        
        % Get to next subplot
        subplot(3,3, (thetaIdx-1) * 3 + flagIdx)
        
        % Create initial discretization data vectors
        % ------------------------------------------
        Xh = linspace(0, 1, n);
        Th = linspace(t0, t, m);
        Uh = Xh' * 0;

        % Get PDE
        [Ah, Mh, Fh] = a07ex01getPDE(         Xh,... domain
                                               f,... rhs
                                           Const,... [a b]
                                           theta,... theta
                                     Th(2)-Th(1),... tau
                                      char(flag)); % FLAG

        % Loop over all time steps
        for t = Th
            % Solve PDE numerical
            Uh = Ah \ (Mh * Uh + (Th(2)-Th(1)) * Fh);

        end % loop over all time steps

        % Solve PDE analyticaly
        U = Xh - (exp(-(1-Xh)/Const(1)) - exp(-1/Const(1))) / (1 - exp(-1/Const(1))); U = U';
        
        % Calculate the maximum error
        err = max(abs(U - Uh));
        
        % Plot both solutions
        % -------------------
        hold on
        plot(Xh,  U, '-r', 'LineWidth', 1.5)
        plot(Xh, Uh, '-k', 'LineWidth', 1.0)
        hold off
        
        % Set axis values
        % ---------------
        grid on
        set(gca, 'FontSize', 15)
        xlabel('X_h', 'FontSize', 12)
        ylabel('u ', 'FontSize', 12)
        legend('Ana', 'Num', 'location', 'southeast')
        title(['Theta = ' num2str(theta) ', Flag = ' char(flag) ', err = ' num2str(err, '%4.1d')], 'FontSize', 15)   
        
        drawnow
        
        % Stability measure ( must be greater than 1, II.64)
        stab = 1 - 2 * (1 - theta) * (Th(2)-Th(1)) / (Xh(2)-Xh(1))^2;
        
        % Calculate and print maximum norm
        fprintf('Theta: %d, m: %i, Stability = %d, Error: %d\n', theta, m, stab, err);
        
    end % of loop over all difference stencils
    
end % of loop over all theta schemes

% Remark
disp('THE REULTS WERE DISCUSSED IN THE ASSINGMENT PDF!')

% Save figure
%export_fig('../Documentation/Figures/a07ex01f.png')