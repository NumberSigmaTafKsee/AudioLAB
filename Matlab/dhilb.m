% dhilb.m - ideal lowpass Hilbert transformer FIR filter
%
% h = dhilb(wc, N) = row vector
%
% N = 2M+1 = filter length (odd)
% wc = cutoff frequency in [rads/sample]

function h = dhilb(wc, N)

M = (N-1)/2;

for k = -M:M,
    if k == 0,
        h(k+M+1) = 0;
    else
        h(k+M+1) = (1 - cos(wc * k)) / (pi * k);
    end
end
