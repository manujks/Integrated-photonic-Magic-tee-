% funinvert - Inverts monotonic function in domain and range intervals (0,1)
%           - Interpolates to find points not specified
%
% Syntax:   xdesired = funinvert(x, y, ydesired)
%
% Inputs:   x,y         - equal length vectors representing sampled function
%           ydesired    - arbitrary vector of values for which x is desired
% Output:   xdesired    - vector of desired x values corresponding to given ydesired
%
% Mar 30, 2005 - MP

function xdes = funinvert(x, y, ydes)

warning off MATLAB:fzero:UndeterminedSyntax
for k = 1:length(ydes)
    if (ydes(k) >= 1)
        xdes(k) = 1;
    elseif (ydes(k) <= 0)
        xdes(k) = 0;
    else
        xdes(k) = fzero(inline('interp1(x, y, xdes) - ydes','xdes','x','y','ydes'),[0 max(x)],[],x,y,ydes(k));
    end
end
warning on MATLAB:fzero:UndeterminedSyntax
