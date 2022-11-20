% Generate 3D taper structure - tapers ONLY the small index-block-matrices
% and their dimensions (unlike in 2D case where the actual index
% distribution is generated)
% MP, Mar 21, 2005
%
% Syntax:   [nlyrso, dlyrsxo, dlyrsyo, zo] = gentaper3d(nlyrs, dlyrsx, dlyrsy, z, interpfcns)
% MxNxP structure in (y,x,z)
%
% Input:    nlyrs  - {P}x[MxN] cell array of P>=2 matrices of the indices of
%                              regions of the cross-section
%           dlyrsx - {P}x[N]   cell array of layer widths for P cross-sections
%                              of N x-layers each (min P>=2 = start and end)
%           dlyrsy - {P}x[M]   cell array of P sets of M y-layers
%           z      - {P-1 sets} points at which to find index distribution
%           interpfcns [P-1] - interpolation (taper) function(s): (0,1) -> (0,1)
%               e.g. @sin, 'milosfcn', inline('sin(x)'),...
%               DEFAULT: should be linear..
%           zslice  - [OPTIONAL] z-coordinate at which to output index distrib. slice
%
% Output:   nlyrso           - [MxNxP] set of interpolated cross-sections
%           dlyrsxo, dlyrsyo - [NxP]   set of interpolated block widths/heights at
%                                      given cross-sections
%           z                - z-coordinate of given cross-sections
%           nslice           - [OPTIONAL] index distribution of 2D slice of structure
%
% TO DO LIST:
% - add auto-determination of # sections along z (to have 1 pixel steps);
% then user can specify more dense (x2, x3...) where the index averaging
% makes different "oversampled" layers different
% (relevant if "to-be-used-later" discretization info can be passed in here)

% Code updates:
% =============
% Mar 21, 2005 - First writing.
% Dec 12, 2006 - Updated help text only.


function [nlyrso, dlyrsxo, dlyrsyo, z, nslice] = gentaper3d(nlyrs, dlyrsx, dlyrsy, z, interpfcns, zslice)
%MM = length(dlyrsx); NN = length(dlyrsy);
if (~iscell(dlyrsx) && ~iscell(dlyrsy))  error('dlyrsx and/*or* dlyrsy must be cell arrays of P>=2 sets of layer thickness vectors.');  end
P = max(length(dlyrsx),length(dlyrsy));
if (P<2)  error('Must have at least two sets (start,end) in one of dlyrsx, dlyrsy.');  end
if (P>=2 && length(dlyrsx) > 1 && length(dlyrsy) > 1 && length(dlyrsx) ~= length(dlyrsy)) ...
        error('dlyrsx and dlyrsy must have equal number of cells or one must have only one.');  end
if (length(dlyrsx) == 1)  for k = 2:P, dlyrsx{k} = dlyrsx{1}; end;  end     % # Layers
if (length(dlyrsy) == 1)  for k = 2:P, dlyrsy{k} = dlyrsy{1}; end;  end


for ksec = 1:P-1
    if iscell(z)
        znorm = (z{ksec} - min(z{ksec}))/(max(z{ksec}) - min(z{ksec}));     % Normalized z-coordinate from cell array
    else    % assume equal length horizontal z-vectors for each section stacked vertically in an array, z.
        znorm = (z(ksec,:) - min(z(ksec,:)))/(max(z(ksec,:)) - min(z(ksec,:)));     % Normalized z-coordinate from matrix
    end

    fcn = interpfcns{ksec};
    dlyrsxo = []; dlyrsyo = []; nlyrso = [];
    znorm = (znorm(1:end-1) + znorm(2:end)) / 2;   % [MP] Find index in center of sections, half way between z values of each two interfaces
    for nz = 1:length(znorm)
        fint = feval(fcn, znorm(nz));     % Interpolation fcn value
        dlyrsxo(:,end+1) = (dlyrsx{ksec} * (1-fint) + dlyrsx{ksec+1} * (fint)).'; % Interpolated widths
        dlyrsyo(:,end+1) = (dlyrsy{ksec} * (1-fint) + dlyrsy{ksec+1} * (fint)).'; % Interpolated heights
        if (iscell(nlyrs) & length(nlyrs) > 1)   nlyrso(:,:,end+1) = (nlyrs{ksec}  * (1-fint) + nlyrs{ksec+1}  * (fint));  end % Interpolated indices
    end
    if (iscell(nlyrs) & length(nlyrs) > 1)   nlyrso = nlyrso(:,:,2:end);  end     % Interpolated indices
end
if (~iscell(nlyrs) || (iscell(nlyrs) && length(nlyrs) == 1))  nlyrso = nlyrs;  end       % [MP] TO DO: make sure also to convert a single-cell cellarray to a normal vector (doesn't happen here)

if iscell(z)
    zz = [];
%    for k = 1:length(z),  zz = [zz z{k}];  end
    for k = 1:length(z),  zz = [zz (z{k}(1:end-1) + z{k}(2:end))/2];  end     % Leave out last overlap point! for case when znorm is at center of segment
    z = zz(:).';
else
%    z = z.'; z = z(:).';         % Need to unwrap z for multiple sections to correspond to joined up index matrix
    z = (z(:,1:end-1) + z(:,2:end))/2; z = z(:).';         % [MP] when znorm is averaged to the segment centers
end
