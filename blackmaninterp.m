% Interpolate function from saved function points
% Used for the warping function for generating Blackman-window coupling
% tapers, etc.
%
% NOTE: Data file must have variables 'x' and 'y', for x and y(x).
%
% Milos Popovic, Mar 24, 2005

function y = f(x,datafile)

S=load(datafile);
y = interp1(S.x, S.y, x, 'spline','extrap');
