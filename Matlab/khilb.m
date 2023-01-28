% khilb.m - lowpass FIR Hilbert trasformer design using Kaiser window.
%
% h = khilb(fs, fc, Df, A)
%
% fc = cutoff frequency in [Hz]
% Df = transition width in [Hz]
% A = stopband ripple attenuation in [dB]
% dhilb(wc, N) = ideal FIR Hilbert transformer

function h = khilb(fs, fc, Df, A)

wc = 2 * pi * fc / fs;
DF = Df / fs;

[alpha, N] = kparm(DF, A);
h = dhilb(wc, N) .* kwind(alpha, N);
