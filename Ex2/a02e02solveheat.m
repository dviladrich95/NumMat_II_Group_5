function U = a02e02solveheat(mode, XPOINT, t, a, c)
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
%
%                                                                 Solution
% -------------------------------------------------------------------------
% Create n array
N = [1 : 50]';

% Boundary condition modes
switch mode
    % Condition for dirichlet boundary conditions
    case 1
        % Calculate the B vector
        B = 16./ (pi * (4 * N - N.^3));

        % Adjust B vector for even indices
        % --------------------------------
        EvenIdx    = rem(N, 2) == 0;
        B(EvenIdx) = 0;
        
        % Loop over all x values
        for xIdx = 1 : length(XPOINT)
            % Calculate the partial solutions
            % -------------------------------
            X = B .* sin(N * pi * XPOINT(xIdx) / sqrt(a));
            T = exp(-(N.^2 * pi^2 + c) * t);
            
            % Solve the heat equation with respect to dirichlet boundary conditions
            U(xIdx, 1) = sum(X .* T);
        
        end % of loop over all x values
        
    % Condition for neuman boundary conditions
    case 2
        % Loop over all x values
        for xIdx = 1 : length(XPOINT)
            % The solution takes place at N = 2
            N = 2;
            
            % Calculate the partial solutions
            % -------------------------------
            X = 1 * cos(N * pi * XPOINT(xIdx) / sqrt(a));
            T = exp(-(N.^2 * pi^2 + c) * t) * exp(-c * t);
            
            % Solve the heat equation with respect to dirichlet boundary conditions
            U(xIdx, 1) = X .* T;
            
        end % of loop over all x values
        
    % Invalid mode input
    otherwise
        % Log error
        errordlg('Invalid mode input!', 'a02e02solveheat', 'modal');
        
end % of boundary condition modes