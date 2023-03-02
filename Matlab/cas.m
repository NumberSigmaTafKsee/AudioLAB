% cas.m - filtering by cascade of 2nd order sections
%
% [y, W] = cas(K, B, A, W, x)
%
% B = Kx3 numerator matrix
% A = Kx3 denominator matrix
% W = Kx3 state matrix
% x = scalar input
% y = scalar output
% based on cas.c

function [y, W] = cas(K, B, A, W, x)

y = x;
       
for i = 1:K,
       [y, W(i,:)] = sos(B(i,:), A(i,:), W(i,:), y);
end
