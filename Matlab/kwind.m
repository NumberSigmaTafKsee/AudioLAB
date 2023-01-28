% kwind.m - Kaiser window.
%
% w  = kwind(alpha, N) = row vector
%
% alpha = Kaiser window shape parameter
% N  = 2M+1 = window length (must be odd)

function w = kwind(alpha, N)

M = (N-1) / 2; 
den = I0(alpha);

for n = 0:N-1,
       w(n+1) = I0(alpha * sqrt(n * (N - 1 - n)) / M) / den;
end
