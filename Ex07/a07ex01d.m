function [Uh, err] = a07ex01d(N, M, T, A_H, M_H, F_H)
% ------------------------------------------------------------------------\
% Assignment 67, Exercise 1d                                              |
%                                                             submitted by|
%                                                                         |
%                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
%                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
%                   Riccardo Parise | 412524 | Scientific Computing, M.Sc.|
%        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
%                                                                         |
%                                                        in  MATLAB R2014a|
% ------------------------------------------------------------------------/
%
%                                                                 Solution
% -------------------------------------------------------------------------
% Generate 1D arrays
% ------------------
Xh = linspace(0, 1, N);
Th = linspace(1E-5, T, M);

% Calculate timestep
tau = Th(2) - Th(1);

% Create intial state for the Heaviside function
Uh(Xh > 1/2) = 1; Uh = Uh';

% Loop over all time steps
for t = Th
    % Solve PDE numerical
    Uh = A_H \ (M_H * Uh + tau * F_H);
        
end % of loop over all time steps

% Solve PDE analytical
U = a07ex01getsol(Xh, t);

% Calculate the error using maximum norm
err = max(abs(U - Uh));

%                                                                     Plot
% -------------------------------------------------------------------------
% Plot results
% ------------
hold on
plot(Xh,  U, '-r','LineWidth', 1.5)
plot(Xh, Uh, '-k','LineWidth', 1.0)
hold off

% Set axis values
% ---------------
grid on
set(gca, 'FontSize', 15)
xlabel('X_h', 'FontSize', 12)
ylabel('u ', 'FontSize', 12)
legend('Ana', 'Num', 'location', 'southeast')
title(['n = ' num2str(N) ', m = ' num2str(M) ', err = ' num2str(err, '%4.1d')], 'FontSize', 15)    

drawnow