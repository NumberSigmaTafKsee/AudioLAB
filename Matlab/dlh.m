% dlh.m - ideal lowpass/highpass FIR filter
%
% h = dlh(s, wc, N) = row vector
%
% s = 1, -1 = lowpass, highpass
% N = 2M+1 = filter length (odd)
% wc = cutoff frequency in [rads/sample]

function h = dlh(s, wc, N)

M = (N-1)/2;

for k = -M:M,
    if k == 0,
        h(k+M+1) = (1-s) / 2 + s * wc / pi;
    else
        h(k+M+1) = s * sin(wc * k) / (pi * k);
    end
end
