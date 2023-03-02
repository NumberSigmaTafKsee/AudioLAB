% kparm2.m - Kaiser window parameters for spectral analysis.
%
% [alpha, L] = kparm2(DF, R)
%
% alpha = window shape parameter
% L = window length (odd)
% DF = Df/fs = mainlobe width in units of fs
% R = relative sidelobe level in dB
% R must be less than 120 dB.

function [alpha, L] = kparm2(DF, R)

c = 6 * (R + 12) / 155;

if R < 13.26
       alpha = 0;
elseif R < 60
       alpha = 0.76609 * (R - 13.26)^0.4 + 0.09834 * (R - 13.26);
else
       alpha = 0.12438 * (R + 6.3);
end

L = 1 + ceil(c / DF);
L = L + 1 - rem(L, 2);                           % next odd integer
