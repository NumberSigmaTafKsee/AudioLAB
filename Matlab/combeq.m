% combeq.m - periodic comb/notch EQ filter design
%
% [a, b, c, beta] = combeq(G0, G, GB, D, Dw, s)
%
% s = 1, -1 for peaks at (2*k)*pi/D, or, (2*k+1)*pi/D
% G0, G, GB = reference, boost/cut, and bandwidth gains
% D = period, Dw = width in [rads/sample]
% note Dw < pi/D
%
% for plain COMB use:  G0=0, G=1, GB=1/sqrt(2)
% for plain NOTCH use: G0=1, G=0, GB=1/sqrt(2)
%
% H(z) = (b - c z^(-D)) / (1 - a z^(-D))   (caution: note minus signs)

function [a, b, c, beta] = combeq(G0, G, GB, D, Dw, s)

beta = tan(D * Dw / 4) * sqrt(abs((GB^2 - G0^2) / (G^2 - GB^2)));
a = s * (1 - beta) / (1 + beta);
b = (G0 + G * beta) / (1 + beta);
c = s * (G0 - G * beta) / (1 + beta);
