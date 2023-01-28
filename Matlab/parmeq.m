% parmeq.m - second-order parametric EQ filter design
%
% [b, a, beta] = parmeq(G0, G, GB, w0, Dw)
%
% b = [b0, b1, b2] = numerator coefficients
% a = [1,  a1, a2] = denominator coefficients
% G0, G, GB = reference, boost/cut, and bandwidth gains
% w0, Dw = center frequency and bandwidth in [rads/sample]
% beta = design parameter
%
% for plain PEAK use:  G0=0, G=1, GB=1/sqrt(2)
% for plain NOTCH use: G0=1, G=0, GB=1/sqrt(2)

function [b, a, beta] = parmeq(G0, G, GB, w0, Dw)

beta = tan(Dw/2) * sqrt(abs(GB^2 - G0^2)) / sqrt(abs(G^2 - GB^2));
b = [(G0 + G*beta), -2*G0*cos(w0), (G0 - G*beta)] / (1+beta);
a = [1, -2*cos(w0)/(1+beta), (1-beta)/(1+beta)];
