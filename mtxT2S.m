%---------.---------.---------.---------.---------.---------.---------.--------|
% Code:      Transfer-to-Scattering Matrix Converter
% Date:      July 2, 2002
% Author:    Milos Popovic
%
% [S(1:2,1:2,:)] = mtxT2S( T(1:2,1:2,:) )
%
% A tensor containing any number of 2x2 T (transfer) matrices is converted
% to a tensor of S (scattering) matrices.
%
% The top and bottom of the output vector of the transfer matrix become the
% input and output for the top row of the scattering matrix.

% Jul  2, 2002 - First written.
%---------.---------.---------.---------.---------.---------.---------.--------|

function [S] = mtxT2S(T)

S(2,1,:) = 1./T(1,1,:);
S(2,2,:) = -T(1,2,:)./T(1,1,:);
S(1,1,:) = T(2,1,:)./T(1,1,:);
S(1,2,:) = T(2,2,:) - T(1,2,:).*T(2,1,:)./T(1,1,:);
