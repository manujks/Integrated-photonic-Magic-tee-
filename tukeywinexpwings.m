function w = tukeywinexpwings(n,r,rwings,power)
%TUKEYWIN Tukey window, with exponential wings at ends (optionally weighted by a power.
%   W = TUKEYWIN(N,R) returns an N-point Tukey window in a column vector.
%   A Tukey window is also known as the cosine-tapered window.  The R 
%   parameter specifies the ratio of taper to constant sections. This ratio 
%   is normalized to 1 (i.e., 0 < R < 1).  Note, at extreme values of R, the 
%   Tukey window degenerates into other common windows. Thus when R = 1, 
%   it is equivalent to a Hanning window. Conversely, for R = 0 the Tukey 
%   window assumes a constant value (i.e., boxcar).
%
%   EXAMPLE:
%      N = 64; 
%      w = tukeywin(N,0.5); 
%      plot(w); title('64-point Tukey window, Ratio = 0.5');
%
%   See also BARTLETT, BARTHANNWIN, BLACKMAN, BLACKMANHARRIS, BOHMANWIN, 
%            CHEBWIN, GAUSSWIN, HAMMING, HANN, KAISER, NUTTALLWIN, RECTWIN,
%            TRIANG, WINDOW.

%   Reference:
%     [1] fredric j. harris [sic], On the Use of Windows for Harmonic Analysis
%         with the Discrete Fourier Transform, Proceedings of the IEEE,
%         Vol. 66, No. 1, January 1978, Page 67, Equation 38.

%   Author(s): A. Dowd
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.9 $  $Date: 2002/03/28 17:31:09 $

error(nargchk(2,4,nargin));

if(nargin < 3) rwings = 0.2; end      % Default exponential wings = 10% on each side
if(nargin < 4) power = 2; end      % Default taper wings up to Gaussian (pwr 2).. to keep exponential (no tapering), set to 1

if r <= 0,
    w = ones(n,1);
elseif r >= 1,
    w = hann(n);
else
    t = linspace(0,1,n)';
    % Defines period of the taper as 1/2 period of a sine wave.
    per = r/2; 
    tl = floor(per*(n-1))+1;
    th = n-tl+1;
    % Window is defined in three sections: taper, constant, taper
%    w = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end) - 1 + per)))/2)];
    w = [ sin(pi/2 * t(1:tl)/per).^2;  ones(th-tl-1,1); flipud(sin(pi/2 * t(1:tl)/per).^2) ];

    % Put in exponential wings with matched slope at contact point
    fracexp = rwings/2;
    tl = floor(fracexp*(n-1))+1;
    th = n-tl+1;
    slope = pi/per * sin(pi/2 * t(tl)/per) * cos(pi/2 * t(tl)/per); % Slope at join point
    
    slopevec = slope * (1 + (power-1) * (1-[0:(tl-1)].'/(tl-1)));
    w(1:tl) = w(tl)*exp(slopevec./w(tl).*(t(1:tl)-t(tl)));
    w(th:end) = flipud(w(1:tl));

end


% [EOF]
