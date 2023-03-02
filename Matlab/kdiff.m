% kdiff.m - lowpass FIR differentiator design using Kaiser window.
%
% h = kdiff(fs, fc, Df, A)
%
% fc = cutoff frequency in [Hz]
% Df = transition width in [Hz]
% A = stopband ripple attenuation in [dB]
% ddiff(wc, N) = ideal FIR differentiator

function h = kdiff(fs, fc, Df, A)

wc = 2 * pi * fc / fs;
DF = Df / fs;

[alpha, N] = kparm(DF, A);
h = ddiff(wc, N) .* kwind(alpha, N);
