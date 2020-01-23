% Sample solution DP
function [ne,np,e1]=generateelements1D(x)

np = length(x);         % number of points
ne = np-1;              % number of elements
e1 = [1:(np-1);2:np]';  % element specificiation