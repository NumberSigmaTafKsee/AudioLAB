% cfir2.m - FIR filter using circular delay-line buffer
%
% [y, w, q] = cfir2(M, h, w, q, x)
%
% h = order-M filter (row vector)
% w = circular filter state (row vector)
% q = circular index into w
% x = scalar input
% y = scalar output
% based on cfir2.c

function [y, w, q] = cfir2(M, h, w, q, x)

w(q+1) = x;                                 % read input
y = h * w(rem(q+(0:M), M+1)+1)';            % compute output
q = cdelay2(M, q);                          % update delay
