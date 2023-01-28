% dtft.m - DTFT of a signal at a frequency vector w
%
% X = dtft(x, w);
%
% x = row vector of time samples
% w = row vector of frequencies in rads/sample
% X = row vector of DTFT values
%
% based on and replaces both dtft.c and dtftr.c

function X = dtft(x, w)

[L1, L] = size(x);

z = exp(-j*w);

X = 0;
for n = L-1:-1:0,
       X = x(n+1) + z .* X;
end
