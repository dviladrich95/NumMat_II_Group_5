function V = a02e01solution(X, Y)
% ------------------------------------------------------------------------\
% Assignment 2, Exercise 1                                                |
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
% Calculate the radius
R = sqrt(X.^2 + Y.^2);

% Adjust the edge data
% --------------------
R (:,  1) = 1;
R (:,end) = 2;

% Condition for all radius elements within the range of [1,2]
if isempty(find(R<1, 1)) && isempty(find(2<R, 1))
    % Calculate the angle
    PHI = atand(Y./X);
    
    % Calculate solution of the PDE
    V = R .* cosd(PHI);

% Exceptional conditions
else
    % Log error exception
    errordlg('The radius is off range. Check X and Y data!', 'a02e01solution', 'modal');

end % of data evaluation conditions