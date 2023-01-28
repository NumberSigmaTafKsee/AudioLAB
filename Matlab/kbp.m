% kbp.m - bandpass FIR filter design using Kaiser window.
%
% h = kbp(fs, fpa, fpb, fsa, fsb, Apass, Astop, s)
%
% s = 1, -1 = standard, alternative design
% dbp(wa, wb, N) = ideal bandpass FIR filter

function h = kbp(fs, fpa, fpb, fsa, fsb, Apass, Astop, s)

Df = min(fpa-fsa, fsb-fpb);  DF = Df / fs;
fa = ((1+s) * fpa + (1-s) * fsa - s * Df) / 2;  wa = 2 * pi * fa / fs;
fb = ((1+s) * fpb + (1-s) * fsb + s * Df) / 2;  wb = 2 * pi * fb / fs;

dpass = (10^(Apass/20) - 1) / (10^(Apass/20) + 1);
dstop = 10^(-Astop/20);
d = min(dpass, dstop);
A = -20 * log10(d);

[alpha, N] = kparm(DF, A);
h = dbp(wa, wb, N) .* kwind(alpha, N);
