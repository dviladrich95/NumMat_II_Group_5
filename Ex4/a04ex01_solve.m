% Workspace commands
clc
clear
close all

% Set domain parameters
N = 1022;
L = 1;

% Get discretization data
[Lh, Xh, Yh, Idx_Bd] = a04ex01_Lh5(L, N);

% Determine the right hand side
RHS = sin(pi*Xh(:)).*sin(pi*Yh(:));

% Determine numerical U 
Uh = Lh\RHS;

% Determine analytical U
U = RHS ./ (2*pi^2);

% Calcuate the maximum norm
errMax = max(max(abs(Uh - U)));

% Plot settings
% -------------
surf(Xh, Yh, reshape(Uh, [N + 2 N + 2]),'EdgeAlpha', 0.1)
colorbar
figScaleFac = 0.75;
grid on
xlabel('x_h', 'FontSize', 15)
ylabel('y_h', 'FontSize', 15)
zlabel('u_h', 'FontSize', 15)
title(['Number of DOFs: ' num2str(length(Lh)) ', Maximum Error: ' num2str(errMax)], 'FontSize', 15)
set(gca, 'FontSize', 15)
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
     
export_fig('../Documentation/Figures/a04ex01Laplace.png')