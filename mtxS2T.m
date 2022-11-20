%---------.---------.---------.---------.---------.---------.---------.--------|
% Code:      Scattering-to-Transfer Matrix Converter
% Date:      July 2, 2002
% Author:    Milos Popovic
%
% [T(1:2,1:2,:)] = mtxS2T( S(1:2,1:2,:) )
%
% A tensor containing any number of 2x2 S (scattering) matrices is converted
% to a tensor of T (transfer) matrices.
%
% The top row input and output for the scattering matrix become the top and
% bottom output for the transfer matrix.

% Jul  2, 2002 - First written.
%---------.---------.---------.---------.---------.---------.---------.--------|

function [T] = mtxS2T(S)

T(1,1,:) = 1./S(2,1,:);
T(1,2,:) = -S(2,2,:)./S(2,1,:);
T(2,1,:) = S(1,1,:)./S(2,1,:);
T(2,2,:) = S(1,2,:) - S(2,2,:).*S(1,1,:)./S(2,1,:);
