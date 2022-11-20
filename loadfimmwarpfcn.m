% Load [x, y(x)] vectors, each from 0..1, from a file formatted as a FIMMWAVE-format
% warp function.
%
% Syntax:   [x,y] = loadfimmwarpfcn(filename)
%
% Milos Popovic, Apr 17, 2005

function [x,y] = loadfimmwarpfcn(filename)
global VMODE; if isempty(VMODE) VMODE = 1; end          % Default chattiness mode

fid = fopen(filename, 'r');
if (fid < 0)
    error(['Can''t open file ', filename]);
end

fname = fscanf(fid, '%s\n',1);
if(VMODE > 1)  fprintf(['Reading data set:  ', fname, '\n\n']);  end
xy = fscanf(fid, '%f %f\n', [2 inf]).';
x = xy(:,1); y = xy(:,2);

fclose(fid);
