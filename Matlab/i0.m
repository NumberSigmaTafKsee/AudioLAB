% I0.m - modified Bessel function of 1st kind and 0th order.
%
% S = I0(x)
%
% defined only for scalar x >= 0
% based on I0.c

function S = I0(x)

eps = 10^(-9);
n = 1; S = 1; D = 1;

while D > (eps * S),
        T = x / (2*n);
        n = n+1;
        D = D * T^2;
        S = S + D;
end
