% lhcheb1.m - lowpass/highpass Chebyshev type 1 filter design
%
% [A, B, P] = lhcheb1(s, fs, fpass, fstop, Apass, Astop)
%
% s = 1, -1 = lowpass, highpass
% design parameters: 
% P = [Wpass, Wstop, epass, estop, Nex, N, f3, a];
% A, B are Kx3 matrices of cascaded second-order sections

function [A, B, P] = lhcheb1(s, fs, fpass, fstop, Apass, Astop)

Wpass = tan(pi * fpass / fs);  Wpass = Wpass^s;
Wstop = tan(pi * fstop / fs);  Wstop = Wstop^s;

epass = sqrt(10^(Apass/10) - 1);
estop = sqrt(10^(Astop/10) - 1);

Nex = acosh(estop/epass) / acosh(Wstop/Wpass);
N = ceil(Nex);  r = rem(N,2);  K = (N - r) / 2;

a = asinh(1/epass) / N;
W3 = Wpass * cosh(acosh(1/epass)/N);
f3 = (fs/pi) * atan(W3^s);                           % 3dB frequency
P = [Wpass, Wstop, epass, estop, Nex, N, f3, a];
W0 = sinh(a) * Wpass;

if r==1,
    G = W0 / (1 + W0);
    A(1,:) = [1, s*(2*G-1), 0];
    B(1,:) = G * [1, s, 0];
else
    G = 1 / sqrt(1 + epass^2);
    A(1,:) = [1, 0, 0];
    B(1,:) = G * [1, 0, 0];
end

for i=1:K,
    th = pi * (N - 1 + 2 * i) / (2 * N);
    Wi = Wpass * sin(th);
    D = 1 - 2 * W0 * cos(th) + W0^2 + Wi^2;
    G = (W0^2 + Wi^2) / D;
    a1 = 2 * (W0^2 + Wi^2 - 1) / D;
    a2 = (1 + 2 * W0 * cos(th) + W0^2 + Wi^2) / D;
    A(i+1,:) = [1, s*a1, a2];
    B(i+1,:) = G * [1, s*2, 1];
end
