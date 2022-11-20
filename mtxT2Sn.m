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
% Mar 18, 2005 - Generalized from 2x2(xM) to 2Nx2N(xM) matrices (T and S)
%---------.---------.---------.---------.---------.---------.---------.--------|

function [S] = mtxT2S(T)
[M,N,P] = size(T); M=M/2; N=N/2;
if (round(M) ~= M || M ~= N)  error('MTXS2T: Input must be a 2Nx2N matrix.');  end

for pp = 1:P
    iT11 = inv(T(1:N,1:N,pp));

    S(N+[1:N],1:N,pp)     = iT11;                                       % S21
    S(N+[1:N],N+[1:N],pp) = -iT11 * T(1:N,N+[1:N],pp);                  % S22
    S(1:N,1:N,pp)         = T(N+[1:N],1:N,pp) * iT11;                   % S11
    S(1:N,N+[1:N],pp)     = T(N+[1:N],N+[1:N],pp) - T(N+[1:N],1:N,:) * iT11 * T(1:N,N+[1:N],pp);    % S12
end
