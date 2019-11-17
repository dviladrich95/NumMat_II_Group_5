function a04ex01_solve(FLAG)
% ------------------------------------------------------------------------\
% Assignment 4, Exercise 1d & 1e                                          |
%                                                             submitted by|
%                                                                         |
%                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
%                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
%                   Riccardo Parise | 412524 | Scientific Computing, M.Sc.|
%        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
%                                                                         |
%                                                        in  MATLAB R2014a|
% ------------------------------------------------------------------------/
% READ FIRST:
% To run the code, please set values as stated below for 'FLAG'
%   * FULL: Full matrix,
%   * RED : reduced matrix

close all
%                                                                 Solution
% ------------------------------------------------------------------------

try
% Set domain parameters
% ---------------------
N = 510;
L = 1;

% Get discretization data
[Lh, Xh, Yh, Idx_Bd] = a04ex01_Lh5(L, N);

% Determine the right hand side
RHS = sin(pi*Xh(:)+pi/8).*sin(pi*Yh(:));

% Start time count 
tic;

% Differecen matrix conditions
switch FLAG
    %                                                            Full form
    % ---------------------------------------------------------------------
    case 'FULL'
        % Get number of active DOFs
        numberActiveDOF = size(Lh,2);
        
        % Determine numerical U 
        Uh = Lh\RHS;
    %                                                         Reduced form
    % ---------------------------------------------------------------------
    case 'RED'       
        % Generate master data vector with number of DOFs
        Uh = zeros(size(Lh,1),1);
                
        % Apply Dirichlet boundary conditions
        Uh(Idx_Bd, 1) = RHS(Idx_Bd);
        
        % Modify the right hand side with boundary conditions
        RHS_mod = RHS - Lh * Uh;
        
        % Determine the active dofs other than boundary points
        ActiveDOF = setdiff((1:length(Lh)), Idx_Bd)';

        % Get number of active DOFs
        numberActiveDOF = length(ActiveDOF);
        
        % Calculate the active data
        Uactive = Lh(ActiveDOF, ActiveDOF)\RHS_mod(ActiveDOF);

        % Assign the active data to the master data vector
        Uh(ActiveDOF,1) = Uactive;
        
    %                                                            Exception
    % ---------------------------------------------------------------------
    otherwise
        ME = MException('a04ex01d:InvalidFlag', 'Invalid diffFlag input!');
        throw(ME)

end % of difference matrix conditions

% End time count
time = toc;

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
title(['Matrix form:' FLAG ', # Active DOFs: ' num2str(numberActiveDOF) ', Maximum Error: ' num2str(errMax) ', Duration: ' num2str(time) ' s'], 'FontSize', 15)
set(gca, 'FontSize', 15)
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
     
%export_fig(['../Documentation/Figures/a04ex01Laplace_' FLAG '.png'])

% Catch exception
catch exception
    % Condition for flag exception
    if strcmp(exception.identifier, 'a04ex01d:InvalidFlag')
        % Log error exception
        errordlg('Invalid flag input! Check diffFlag.', 'a04ex01d', 'modal');
    
    % Condition for standard exception
    else
        rethrow(exception)
        
    end % of exception conditions
    
end % of try-catch block