% cdelay2.m - circular delay
%
% q = cdelay2(D, q)
%
% q = 0, 1, ..., D = circular index pointing to w(q+1)
% based on cdelay2.c

function q = cdelay2(D, q)

q = q - 1;                  % decrement index and wrap mod-(D+1)
q = wrap2(D, q);            % when q=-1, it wraps around to q=D
