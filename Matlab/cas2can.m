% cas2can.m - cascade to canonical
%
% a = cas2can(A)
%
% convolves the rows of A
% A is Kx3 coefficient matrix (or, Kx5 for bandpass filters)
% based on cas2can.c

function a = cas2can(A)

[K, L] = size(A);

a = [1];
for i=1:K,
       a = conv(a, A(i,:));
end
