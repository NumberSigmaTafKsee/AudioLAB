% sos.m - second order section
%
% [y, w] = sos(b, a, w, x)
%
% b = [b0, b1, b2] = 3-dim numerator
% a = [1,  a1, a2] = 3-dim denominator
% w = 3-dim filter state (row vector)
% x = scalar input
% y = scalar output
% based on sos.c

function [y, w] = sos(b, a, w, x)

w(1) = x - a(2) * w(2) - a(3) * w(3);
y = b(1) * w(1) + b(2) * w(2) + b(3) * w(3);
w(3) = w(2);
w(2) = w(1);
