% Sech window for optical structures - like Tukey window but uses raised
% exponential tails instead of raised cosine
% Milos Popovic, Apr 20, 2005

function w = sechwin(npoints,ratio,sechlength)
if (nargin < 3) sechlength = 5; disp(['sechlength auto-set to: ' num2str(sechlength)]); end

if ratio <= 0,
    w = ones(npoints,1);          % Uniform window (ratio = 0)
else
    t = linspace(0,1,npoints)';
    % Defines period of the taper as 1/2 period of a sine wave.
    sechfrac = ratio/2;
    tl = floor(sechfrac*(npoints-1))+1;
    th = npoints-tl+1;
    % Window is defined in three sections: taper, constant, taper
    w = [ sech(sechlength * (1 - t(1:tl)/t(tl)));  ones(th-tl-1,1); sech(sechlength * (t(th:end)-t(th))/(t(end)-t(th)))];
end