% kparm.m - Kaiser window parameters for filter design.
%
% [alpha, N] = kparm(DF, A)
%
% alpha = window shape parameter
% N = window length (odd)
% DF = Df/fs = transition width in units of fs
% A = ripple attenuation in dB; ripple = 10^(-A/20)

function [alpha, N] = kparm(DF, A)

if A > 21,                                       % compute D factor
       D = (A - 7.95) / 14.36;
else
       D = 0.922;
end

if A <= 21,                                      % compute shape parameter
       alpha = 0;
elseif A < 50
       alpha = 0.5842 * (A - 21)^0.4 + 0.07886 * (A - 21);
else
       alpha = 0.1102 * (A - 8.7);
end

N = 1 + ceil(D / DF);                            % compute window length
N = N + 1 - rem(N, 2);                           % next odd integer
