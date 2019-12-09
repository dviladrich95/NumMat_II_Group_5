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
n =  2^10;
l = 4.1;

% Set domain matrices
% -------------------
Uh  = zeros(n,n);
U   = zeros(n,n);
Rhs = zeros(n,n);
RR = zeros(n,n);
PHIPHI = zeros(n,n);

% Function definitions for domain and boundary conditions
% -------------------------------------------------------
is_in_domain = @(xx,yy,L) (xx-L/2).^2 + (yy-L/2).^2 <= 4 &...
               (xx-L/2).^2 + (yy-L/2).^2 >= 1;
isOuter = @(xx,yy,L) (xx-L/2).^2 + (yy-L/2).^2 > 1.5^2;
isInner = @(xx,yy,L) (xx-L/2).^2 + (yy-L/2).^2 < 1.5^2;
f            = @(xx,yy) 1;
gIn          = @(xx,yy) 0;
gOut         =@(xx,yy) xx;

% Get numerical differentiation data and domain
[Lh        ,...
 XX        ,...
 YY        ,...
 BarOmega_h,...
 Omega_h   ,...
 GammaOut_h,...
 GammaIn_h ] = a05ex01_get_laplace(l,n,is_in_domain, isOuter, isInner);

% Apply initial domain data
%Uh(:,:)           = g(XX, YY);

% Build right hand side with respect to boundary conditions
% ---------------------------------------------------------
Rhs(Omega_h)   = +f(XX(Omega_h), YY(Omega_h));
Rhs(GammaOut_h)   = +gOut(XX(GammaOut_h)-l/2, YY(GammaOut_h));
Rhs(GammaIn_h)   = +gIn(XX(GammaIn_h)-l/2, YY(GammaIn_h));

% Numerical solution
Uh(BarOmega_h) = Lh \ Rhs(BarOmega_h);

% Analytical solution
RR(BarOmega_h)=sqrt((XX(BarOmega_h)-l/2).^2 + (YY(BarOmega_h)-l/2).^2);
PHIPHI(BarOmega_h)=atan2(YY(BarOmega_h)-l/2,XX(BarOmega_h)-l/2);


U(BarOmega_h) = 2*cos(PHIPHI(BarOmega_h)).*(2/3*(RR(BarOmega_h)-1./RR(BarOmega_h)));




% Calculate the error
errMax = max(max(abs(U(BarOmega_h) - Uh(BarOmega_h))));

% Plot Solution
% -------------
%{-
figure(3)
surf(XX,YY,U-Uh, 'EdgeAlpha', 0.2)
figure(4)
surf(XX,YY,U, 'EdgeAlpha', 0.2)
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
     
%export_fig('../Documentation/Figures/a06ex02.png')
%}