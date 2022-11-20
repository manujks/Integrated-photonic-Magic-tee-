% Film (mode) matching method in 3D - Propagation code
%
% Syntax:   [S] = fmm3d(K,B,z,zscale)
%
% Inputs:   K - (:,:,N-2) set of coupling matrices
%           B - (:,N-1) set of propagation constants
%           z - (N) vector of coordinates for interfaces
%           zscale - (M) vector of any number of scaling values
% Outputs:  S - Scattering matrix set (1:NMODES,1:NMODES,1:length(zscale))
% MP Mar 27, 2005

% Code updates:
% -------------
% Apr 15, 2005 - Added only text description of code inputs above

function [S] = fmm3dprop(K,B,z,zscale)
VMODE = 2;                                                      % Diagnostics (1 = text, 2 = text+plots)
NMODES = size(K,1)/2;
NZ = length(z);

fprintf('\nFMM3Dprop: Z-scaling...\n');
%for mm = 1:length(zscale)
for mm = 1:size(zscale,2)
    T = eye(2*NMODES);                                          % Start with identity matrix
%    zn = z * zscale(mm);
    zn = zscale(:,mm) * z;                                      % Rewritten to permit different optical paths for different modes

    for nz = 1:NZ-1
        sqrtL = spdiags(B(:,nz),0,NMODES,NMODES);               % Propagation constants (root-Eigenvalues) of nz-th section

%        Pnz = [diag(exp(i*diag(sqrtL)*(zn(nz+1)-zn(nz)))) zeros(NMODES); zeros(NMODES) diag(exp(-i*diag(sqrtL)*(zn(nz+1)-zn(nz))))];     % Slice # nz propagator matrix
        Pnz = [diag(exp(i*diag(sqrtL).*(zn(nz+1)-zn(nz)))) zeros(NMODES); zeros(NMODES) diag(exp(-i*diag(sqrtL).*(zn(nz+1)-zn(nz))))];     % Slice # nz propagator matrix
        T = Pnz * T;
        if (nz < NZ-1)
            T = K(:,:,nz) * T;
        end
    end

    T = [T(NMODES+[1:NMODES],NMODES+[1:NMODES]) T(NMODES+[1:NMODES],1:NMODES); T(1:NMODES,NMODES+[1:NMODES]) T(1:NMODES,1:NMODES)]; % Flip around..
    SS = mtxT2Sn(T);
    SS = [SS(NMODES+[1:NMODES],NMODES+[1:NMODES]) SS(NMODES+[1:NMODES],1:NMODES); SS(1:NMODES,NMODES+[1:NMODES]) SS(1:NMODES,1:NMODES)]; % Flip around..
    S(:,:,mm) = SS;
    fprintf('%f.. ', zscale(1,mm)); if(mod(mm,10)==0) fprintf('\n'); end                 % DIAGNOSTICS
    if(VMODE > 1)
        S11(mm) = S(NMODES+1,1,mm); S21(mm) = S(NMODES+1,2,mm);   % Forwards
        S11b(mm) = S(1,NMODES+1,mm); S21b(mm) = S(2,NMODES+1,mm);  % Backwards
        subplot(3,2,5); loglog(zscale(1,1:mm)*max(z), abs([S11(:) S21(:)]).^2, zscale(1,1:mm)*max(z), abs([S11b(:) S21b(:)]).^2, 'o', 'MarkerSize', 3);
        title('S-params vs. length'); axis tight; ylim([1e-6 2]); drawnow;
    end
end
