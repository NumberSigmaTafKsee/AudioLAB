% delay.m - M-fold delay
%
% w = delay(M, w)
%
% w is the (M+1)-dimensional linear delay-line buffer
% based on delay.c

function w = delay(M, w)

w(M+1:-1:2) = w((M+1:-1:2)-1);
