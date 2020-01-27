function Error = a03ex03solve()
% ------------------------------------------------------------------------\
% Assignment 3, Exercise 2c                                               |
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
%                                                                 Solution
% -------------------------------------------------------------------------
Error = [];
H     = [];

% Set array for power value in determining n
P = 1:15;

% Loop over all P
for p = P
    % Get boundary value problem
    [Xh, Lh, Fh] = a03ex03getBVP(p);

    % Calculate max. norm of the error and add to vector
    Uh = Lh \ Fh;
    U = 1 + 4 * Xh.^2 - 3 * Xh.^3;
    Error(end + 1) = max(abs(Uh - U));
                         
    % Calculate step size and add to vector
    H(end + 1) = 1 / length(Xh);
    
end % of loop over all p

% Plot settings
% -------------
plot(H, Error,'LineWidth', 2)
figScaleFac = 0.75;
grid on
xlabel('h'    , 'FontSize', 15)
ylabel('Error', 'FontSize', 15)
set(gca, 'FontSize', 15       ,...
         'XScale'  , 'log'    ,...
         'YScale'  , 'log'    ,...
         'XDir'    , 'reverse')
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
     
%export_fig('../Documentation/Figures/a03ex03Error.png')