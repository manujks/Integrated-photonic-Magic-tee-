% Film (mode) matching method in 3D - (K,B)-Matrix generation code
%
% Syntax:   [K,B,Ex,Ey] = fmm3d(nlyrs, dlyrsx, dlyrsy, dx, dy, k0, {OPTS}, yslice)
%
% Inputs:
%           nlyrs   - [MxNxP] - length P set of MxN index matrices, one for
%                               each slice (generated by gentaper3d.m)
%           dlyrsx  - [NxP]   - length P set of N-vectors of x widths
%           dlyrsy  - [NxP]   - length P set of N-vectors of y heights
%           dx,dy   - discretization
%           k0      - free-space wavelength (2*pi/lambda)
%           OPTS    - [OPTIONAL] OPTS structure for modesolver, can specify:
%                                'mu_guess'    : nco * 2*pi/lambda
%                                'NMODES_CALC' : 2
%                                'tol'         : 1e-8
%           yslice - [OPTIONAL] y-coordinate at which to output 2D
%                                index distrib (slice is at y = dy/2*nyslice)
%
% Outputs:
%           K - Coupling matrices for all interfaces
%           B - Propagation constants of used modes for all sections
% MP Mar 18, 2005

% Code updates:
% -------------
% Aug 30, 2006 - Just modified help text.
% Dec 12, 2006 - Removed some hard-coded values and simplified/generalized
%                input interface to make it more usable.  THIS MODIFICATION
%                BREAKS BACKWARD COMPATIBILITY due to change of input commandline.
% Dec 13, 2006 - Changed nyslice input parameter to yslice (specifying
%                index nyslice was way too cumbersome and hard to use);
%                also now if no yslice is specified as input, but nslice is
%                still given as an output argument, then yslice is assumed
%                to be halfway along the height.

%function [K,B,nslice] = fmm3d(nlyrs, dlyrsx, dlyrsy, dx, dy, k0, NMODES, z, excit, nyslice)   % Pre-20061212 commandline
function [K,B,nslice] = fmm3d(nlyrs, dlyrsx, dlyrsy, dx, dy, k0, OPTS, yslice)    %NMODES, z, excit, nyslice) [20061212]
disp('[fmm3d] Warning: changes made on 20061212 to improve interface may break older scripts (last few input parameters changed).  Please check.');
NZ = size(dlyrsx,2);    % number of sections in z-coordinate
%eta = 4e-7*pi * 299792458;  % Free space impedance = mu * c ~ 377 ohms

%OPTS = struct('mu_guess',1.7*k0,'NMODES_CALC',3,'tol',1e-8);
if (nargin < 7 || isempty(OPTS)) OPTS = struct('mu_guess',max(nlyrs(:))*k0,'NMODES_CALC',3,'tol',1e-8); end   % Default settings
t0 = clock; t1 = clock;
for nz = 1:NZ
    fprintf('Slice %d/%d, slice time = %f sec, elapsed time = %f sec\n', nz, NZ, etime(clock, t1), etime(clock, t0));
    t1 = clock;
    if (size(nlyrs,3) > 1)  nn = nlyrs(:,:,nz);  else  nn = nlyrs(:,:,1);  end;     % Is one index matrix, or one per crosssection specified?
    [N,F,V] = sisolver3d(nn, dlyrsx(:,nz), dlyrsy(:,nz), [dx dy], k0, OPTS);
%    [N] = sisolver3d(nn, dlyrsx(:,nz), dlyrsy(:,nz), [dx dy], k0, OPTS);

    if (nargin < 8) nyslice = round(size(N.n,2)/2); else nyslice = round(yslice/dy/2); end      % Determine where to slice (default: slice halfway along height)
    if(nargout > 2) nslice(:,nz) = N.n(:,nyslice); end              % Store index cross-section at y = yslice (where yslice ~ dy*nyslice)

    % --- Draw computed modes for view ---
    NDRAWMODES_MAX = 9;                             % Maximum # of modes to draw (hard-configurable constant for fmm3d function)
    NDRAWMODES = min( round(length(F.beta)/2 + 0.25) * 2 - 1, NDRAWMODES_MAX );  % Pick all modes (must be odd number), capped at NDRAWMODES_MAX.
    subplot((1+NDRAWMODES)/2,2,1); imagesc(N.n.'); set(gca,'ydir','normal'); %axis image;
    for ip = 1:NDRAWMODES
        subplot((1+NDRAWMODES)/2,2,1+ip); imagesc(F.Ex(:,:,ip).'); set(gca,'ydir','normal'); %axis image;
    end
    drawnow;
    % ------------------------------------
    fprintf('Mode effective indices: '); fprintf('%f, ', F.beta(1:end-1)/k0); fprintf('%f.\n', F.beta(end)/k0);    % Show mode indices

    for k = 1:size(F.Ex,3)  % Power normalize the modes
        P = 1/2 * ecrosshdotz(F,F,[k k],1);                   % Find power
%        exmax = F.Ex(find(abs(F.Ex(:)) == max(max(abs(F.Ex(:,:,k))))));
        cc = 1/sqrt(P); %/(exmax/abs(exmax));
        F.Ex(:,:,k) = F.Ex(:,:,k)*cc; F.Ey(:,:,k) = F.Ey(:,:,k)*cc; F.Ez(:,:,k) = F.Ez(:,:,k)*cc;
        F.Hx(:,:,k) = F.Hx(:,:,k)*cc; F.Hy(:,:,k) = F.Hy(:,:,k)*cc; F.Hz(:,:,k) = F.Hz(:,:,k)*cc;
        fprintf('Power = %f, fields rescaled by cc = %f\n', P, cc);
    end
    B(:,nz) = F.beta;                                           % Save effective index for output

    % Save eigenvectors as guess-vectors for next iteration
%    OPTS.v0 = sum(V(:,:),2);
    OPTS.mu_guess = (F.beta(1));%+F.beta(2))/2;    % For eigs eigensolver (use shift-invert for first computation, then 'LR' with guess vector!)=
%    OPTS.sigma = 'LR';
%    OPTS.v0 = V; OPTS.mu_guess = F.beta;            % For jdqr eigensolver
%    OPTS.mu_guess = max(F.beta(:))*1.05;  % Use max eigenvalue + 5% as guess for next time [20061212]

    if (nz > 1)
        NMODES = size(F.Ex,3);
        for m = 1:NMODES        % Output side modes (next section)
            for k = 1:NMODES    % Input side modes (previous section)
                EkHm = ecrosshdotz(Fold,F,[k m],1); EmHk = conj( ecrosshdotz(F,Fold,[m k],1) );
                K(m,k,nz-1)               = 1/4 * (EkHm + EmHk);    % Save coupling matrices for output (new section modes index dim 1=m, old sec modes index dim 2=k)
                K(m+NMODES,k,nz-1)        = 1/4 * (EkHm - EmHk);
                K(m,k+NMODES,nz-1)        = K(m+NMODES,k,nz-1);  K(m+NMODES,k+NMODES,nz-1) = K(m,k,nz-1);   % [20061212] It is *correct* that these are not transposes of the two
                                                                                                            % matrices on lines above; however, EACH quarter matrix should
                                                                                                            % be symmetric in the limit of small dz discretization, but in
                                                                                                            % reality (coarse dz) they are not since mode patterns differ
                                                                                                            % - this is not ideal.
            end
        end
        %          -- Sign correction of mode fields to keep diagonal coupling positive --
        SG = diag( sign(diag(K(1:end/2,1:end/2,nz-1))) );
        K(:,:,nz-1) = [SG zeros(NMODES,NMODES); zeros(NMODES,NMODES) SG] * K(:,:,nz-1);
        for k = find(diag(SG)<0)
            fprintf('Flipping sign of output facet mode %d...\n',k);
            F.Ex(:,:,k) = -F.Ex(:,:,k);
            F.Ey(:,:,k) = -F.Ey(:,:,k);
            F.Hx(:,:,k) = -F.Hx(:,:,k);
            F.Hy(:,:,k) = -F.Hy(:,:,k);     % Flip transverse fields (ignore longitudinal fields for now.. I think they don't change anyway)
        end
    end
    Fold = F;
end


% OLD COPY OF CODE SNIPPETS TO EVENTUALLY BE REMOVED
% --------------------------------------------------
%
% [Removed on 20061212, line 3:]
% Syntax:   [K,B,Ex,Ey] = fmm3d(nlyrs, dlyrsx, dlyrsy, dx, dy, k0, NMODES, {z, excit, nyslice})
%           nlyrs   - {P}x[MxN] - length P cell array of MxN index
%                     matrices, one for each slice (generated by gentaper3d.m)
%           dlyrsx  - {P}x[1xN] - length P cell array