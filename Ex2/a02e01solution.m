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
Radius = sqrt(X.^2 + Y.^2);

% Adjust the edge data
% --------------------
Radius (:,  1) = 1;
Radius (:,end) = 2;

% Condition for all radius elements within the range of [1,2]
if isempty(find(Radius<1, 1)) && isempty(find(2<Radius, 1))
    % Calculate the angle
    Angle = atan2d(Y,X);
    
    % Calculate the partial terms
    % ---------------------------
    R   = (Radius - 1./Radius);
    PHI = 4/3 * cosd(Angle);
    
    % Calculate solution of the PDE
    V = R .* PHI;

% Exceptional conditions
else
    % Log error exception
    errordlg('The radius is off range. Check X and Y data!', 'a02e01solution', 'modal');

end % of data evaluation conditions