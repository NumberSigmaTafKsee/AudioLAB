% fir.m - sample processing algorithm for FIR filter.
%
% [y, w] = fir(M, h, w, x)
%
% h = order-M filter (row vector)
% w = filter state (row vector)
% x = scalar input
% y = scalar output
% based on fir2.c

function [y, w] = fir(M, h, w, x)

w(1) = x;                          % read input
y = h * w';                        % compute output
w = delay(M, w);                   % update delay
