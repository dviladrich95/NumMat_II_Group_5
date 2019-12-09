% ------------------------------------------------------------------------\
% Assignment 5, Exercise 1d & 1e &1f                                      |
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
%                                                              Solution d)
% -------------------------------------------------------------------------
% Set domain parameters
% ---------------------
n =  2^5;
l = 2.2;

% Set domain matrices
% -------------------
Uh  = zeros(n,n);
U   = zeros(n,n);
Rhs = zeros(n,n);

% Function definitions for domain and boundary conditions
% -------------------------------------------------------
is_in_domain = @(xx,yy,L) (xx-L/2).^2 + (yy-L/2).^2 <= 1;
f            = @(xx,yy) 1;
g            = @(xx,yy) 0;

% Get numerical differentiation data and domain
[Lh        ,...
 XX        ,...
 YY        ,...
 BarOmega_h,...
 Omega_h   ,...
 Gamma_h   ] = a05ex01_get_laplace(l,n,is_in_domain);

% Apply initial domain data
Uh(:,:)           = g(XX, YY);

% Build right hand side with respect to boundary conditions
% ---------------------------------------------------------
Rhs(Omega_h)   = +f(XX(Omega_h), YY(Omega_h));
Rhs(Gamma_h)   = +g(XX(Gamma_h), YY(Gamma_h));

% Numerical solution
Uh(BarOmega_h) = Lh \ Rhs(BarOmega_h);

% Alaytical solution
U(BarOmega_h) = 0.25 -((XX(BarOmega_h)-l/2).^2 + (YY(BarOmega_h)-l/2).^2)/4;

% Calculate the error
errMax = max(max(abs(U(BarOmega_h) - Uh(BarOmega_h))));

% Plot Solution
% -------------
%{-
figure(3)
surf(XX,YY,Uh, 'EdgeAlpha', 0.2)
view([-38, 54])
colorbar
grid on
axis tight equal
figScaleFac = 0.75;
xlabel('x_h', 'FontSize', 15)
ylabel('y_h', 'FontSize', 15)
zlabel('u_h', 'FontSize', 15)
title(['# Active DOFs: ' num2str(length(Lh)) ', Maximum Error: ' num2str(errMax)], 'FontSize', 15)
set(gca, 'FontSize', 15)
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
     
%export_fig(['../Documentation/Figures/a05ex01Laplace_256.png'])
%}
%                                                              Solution e)
% -------------------------------------------------------------------------
%{-
ErrorMax = [];
H        = [];

% Loop over all n
for n = 2.^(2:9)
    % Set domain matrices
    % -------------------
    Uh  = zeros(n,n);
    U   = zeros(n,n);
    Rhs = zeros(n,n);
    
    % Get numerical differentiation data and domain
    [Lh        ,...
     XX        ,...
     YY        ,...
     BarOmega_h,...
     Omega_h   ,...
     Gamma_h   ] = a05ex01_get_laplace(l,n,is_in_domain);

    Rhs(Omega_h)   = +f(XX(Omega_h),YY(Omega_h));
    Rhs(Gamma_h)   = +g(XX(Gamma_h),YY(Gamma_h));

    % Numerical solution
    Uh(BarOmega_h) = Lh \ Rhs(BarOmega_h);

    % Alaytical solution
    U(BarOmega_h) = 0.25 -((XX(BarOmega_h)-l/2).^2 + (YY(BarOmega_h)-l/2).^2)/4;

    % Calculate the error
    ErrorMax(end + 1) = max(max(abs(U(BarOmega_h) - Uh(BarOmega_h))));
    
    % Calculate step size and add to vector
    H(end + 1) = 1 / length(XX);
    
end % of loop over all p

% Calculate the convergence order
order = log10(ErrorMax(end) / ErrorMax(1)) / log10(H(end) / H(1));

% Plot settings
% -------------
figure(4)
plot(H, ErrorMax,'LineWidth', 2)
figScaleFac = 0.75;
grid on
xlabel('h'    , 'FontSize', 15)
ylabel('Error', 'FontSize', 15)
title(['Order of the convergence: ' num2str(order) ]) 
set(gca, 'FontSize', 15       ,...
         'XScale'  , 'log'    ,...
         'YScale'  , 'log'    ,...
         'XDir'    , 'reverse')
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])
     
%export_fig('../Documentation/Figures/a05ex01Error.png')
%}



%                                                              Solution f)
% -------------------------------------------------------------------------
% More complicated domains and problems
% -------------------------------------
isInDomain2 = @(xx,yy,L) (xx - L/12) - 0.05*cos(30 * yy / L)       >= 0 & ...
                         (xx - L/1.1) - 0.1*cos(30 * yy / L)       <= 0 & ...
                         (yy - L/1.09) - 0.1*sin(10 * pi * xx / L) <= 0 & ...
                         (yy - L/12) - 0.1*sin(10 * pi * xx / L)   >= 0 &...
                          0.8*(xx-L/2).^2 + 1.3*(yy-L/2).^2 > 0.1;                     
f2          = @(xx,yy) 5*xx+3*yy;
g2          = @(xx,yy) 0.2*sin(6*pi*xx).*sin(6*pi*yy) + 0.5*sin(2*xx);

% Get numerical differentiation data and domain
[Lh        ,...
 XX        ,...
 YY        ,...
 BarOmega_h,...
 Omega_h   ,...
 Gamma_h   ] = a05ex01_get_laplace(l,n,isInDomain2);

% Apply initial domain data
Uh           = g2(XX/l, YY/l);

% Build right hand side with respect to boundary conditions
% ---------------------------------------------------------
Rhs(Omega_h)   = +f2(XX(Omega_h)  , YY(Omega_h)  );
Rhs(Gamma_h)   = +g2(XX(Gamma_h)/l, YY(Gamma_h)/l);

% Numerical solution
Uh(BarOmega_h) = Lh \ Rhs(BarOmega_h);

% Plot Solution
% -------------
figure(5)
surf(XX,YY,Uh, 'EdgeAlpha', 0.2)
view([-38, 54])
colorbar
grid on
axis tight equal
figScaleFac = 0.75;
xlabel('x_h', 'FontSize', 15)
ylabel('y_h', 'FontSize', 15)
zlabel('u_h', 'FontSize', 15)
set(gca, 'FontSize', 15)
set(gcf, 'Color'            , 'white'                                      ,...
         'PaperSize'        , [34, 34]                                     ,...
         'PaperPositionMode', 'auto'                                       ,...
         'Position'         , [0, 0, 1280 * figScaleFac, 768 * figScaleFac])