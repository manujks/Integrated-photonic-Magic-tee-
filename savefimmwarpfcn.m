% Save [x, y(x)] vectors, each from 0..1, in a file as a FIMMWAVE-format
% warp function.
%
% Syntax:   savefimmwarpfcn(filename, x, y)
%
% Milos Popovic, Apr 14, 2005

function savefimmwarpfcn(filename, x, y)

fid = fopen(filename, 'w');
if (fid < 0)
    error(['Can''t open file ', filename]);
end

fprintf(fid, '%s\n', filename);
for k = 1:length(x)
    fprintf(fid, '%f %f\n', x(k), y(k));
end

fclose(fid);
