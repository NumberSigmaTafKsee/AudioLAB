% wrap2.m - circular wrapping of delay line buffer
%
% q = wrap2(D, q)
%
% q = 0, 1, ..., D = circular index pointing to w(q+1)
% based on wrap2.c

function q = wrap2(D, q)

if q > D,
       q = q - (D+1);              % if q=D+1, it wraps to q=0
end

if q < 0,
       q = q + (D+1);              % if q=-1, it wraps to q=D
end
