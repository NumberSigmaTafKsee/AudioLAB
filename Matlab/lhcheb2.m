% lhcheb2.m - lowpass/highpass Chebyshev type 2 filter design
%
% [A, B, P] = lhcheb2(s, fs, fpass, fstop, Apass, Astop)
%
% s = 1, -1 = lowpass, highpass
% design parameters:
% P = [Wpass, Wstop, epass, estop, Nex, N, f3, a];
% A, B are Kx3 matrices of cascaded second-order sections

function [A, B, P] = lhcheb2(s, fs, fpass, fstop, Apass, Astop)

Wpass = tan(pi * fpass / fs);  Wpass = Wpass^s;
Wstop = tan(pi * fstop / fs);  Wstop = Wstop^s;

epass = sqrt(10^(Apass/10) - 1);
estop = sqrt(10^(Astop/10) - 1);

Nex = acosh(estop/epass) / acosh(Wstop/Wpass);
N = ceil(Nex);  r = rem(N,2);  K = (N - r) / 2;

a = asinh(estop) / N;
W3 = Wstop / cosh(acosh(estop)/N);
f3 = (fs/pi) * atan(W3^s);
P = [Wpass, Wstop, epass, estop, Nex, N, f3, a];
W0 = sinh(a) / Wstop;                              % reciprocal of text

if r==1,
    G = 1 / (1 + W0);
    A(1,:) = [1, s*(2*G-1), 0];
    B(1,:) = G * [1, s, 0];
else
    A(1,:) = [1, 0, 0];
    B(1,:) = [1, 0, 0];
end

for i=1:K,
    th = pi * (N - 1 + 2 * i) / (2 * N);
    Wi = sin(th) / Wstop;                          % reciprocal of text
    D = 1 - 2 * W0 * cos(th) + W0^2 + Wi^2;
    G = (1 + Wi^2) / D;
    b1 = 2 * (1 - Wi^2) / (1 + Wi^2);
    a1 = 2 * (1 - W0^2 - Wi^2) / D;
    a2 = (1 + 2 * W0 * cos(th) + W0^2 + Wi^2) / D;
    A(i+1,:) = [1, s*a1, a2];
    B(i+1,:) = G * [1, s*b1, 1];
end
