% Compute stretch function to get from known coupling distribution k, to
% desired distribution g.
% M. Popovic, Mar 24, 2005
%
% Syntax:   fvec = fstretch(kappavec, gvec)
%
% Inputs:   kappavec    - actual coupling vector, kappa(z)
%           gvec        - desired coupling vector, g(z)
%           z           - OPTIONAL non-uniform z vector (if uniform can be left out)
%
% Output:   fvec        - required coordinate stretching function, f(z)

function fvec = fstretch(kappavec, gvec, z)

gvecAD = cumtrapz(gvec)/trapz(gvec);                            % Antiderivatives
kappavecAD = cumtrapz(kappavec)/trapz(kappavec);
zn = ([1:length(kappavec)]-1)/(length(kappavec)-1); % Normalized coordinate..
%gvecAD = cumtrapz([0; gvec(:); gvec(end)])/trapz([0; gvec(:); gvec(end)]);                            % Antiderivatives
%kappavecAD = cumtrapz([0; kappavec(:); kappavec(end)])/trapz([0; kappavec(:); kappavec(end)]);
%zn = ([1:length(kappavec)+2]-1)/(length(kappavec+2)-1); % Normalized coordinate..

fznguess = 0;
warning off MATLAB:fzero:UndeterminedSyntax
for k = 1:length(zn)
%    fvec(k) = fzero(@zerofcn, fznguess, [], gvecAD(k), zn, kappavecAD, 'spline', 'extrap');
    fvec(k) = fzero(@zerofcn, fznguess, [], gvecAD(k), zn, kappavecAD, 'linear', 'extrap');
    fznguess = fvec(k);
end
warning on MATLAB:fzero:UndeterminedSyntax


% Zeroing function
function y = zerofcn(fznval, gvecADval, znvec, kappavecAD, typeflag, extrapflag)
y = interp1(znvec, kappavecAD, fznval, typeflag, extrapflag) - gvecADval;
