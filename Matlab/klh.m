% klh.m - lowpass/highpass FIR filter design using Kaiser window.
%
% h = klh(s, fs, fpass, fstop, Apass, Astop)
%
% s = 1, -1 = lowpass, highpass
% dlh(s, wc, N) = ideal lowpass/highpass FIR filter

function h = klh(s, fs, fpass, fstop, Apass, Astop)

fc = (fpass + fstop) / 2;  wc = 2 * pi * fc / fs;
Df = s * (fstop - fpass);  DF = Df / fs;

dpass = (10^(Apass/20) - 1) / (10^(Apass/20) + 1);
dstop = 10^(-Astop/20);
d = min(dpass, dstop);  
A = -20 * log10(d);

[alpha, N] = kparm(DF, A);
h = dlh(s, wc, N) .* kwind(alpha, N);
