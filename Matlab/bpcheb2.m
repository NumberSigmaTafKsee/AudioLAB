% bpcheb2.m - bandpass Chebyshev type 2 digital filter design
%
% [A, B, P] = bpcheb2(fs, fpa, fpb, fsa, fsb, Apass, Astop)
%
% design parameters:
% P = [Wpass, Wstop, Wsa, Wsb, c, epass, estop, Nex, N, f3a, f3b, a];
% A, B are Kx5 matrices of cascaded fourth-order sections

function [A, B, P] = bpcheb2(fs, fpa, fpb, fsa, fsb, Apass, Astop)

c = sin(2*pi*(fpa + fpb)/fs)/(sin(2*pi*fpa/fs) + sin(2*pi*fpb/fs));
Wpass = abs((c - cos(2*pi*fpb/fs)) / sin(2*pi*fpb/fs));
Wsa = (c - cos(2*pi*fsa/fs)) / sin(2*pi*fsa/fs);
Wsb = (c - cos(2*pi*fsb/fs)) / sin(2*pi*fsb/fs);
Wstop = min(abs(Wsa), abs(Wsb));

epass = sqrt(10^(Apass/10) - 1);
estop = sqrt(10^(Astop/10) - 1);

Nex = acosh(estop/epass) / acosh(Wstop/Wpass);
N = ceil(Nex);  r = rem(N,2);  K = (N - r) / 2;

a = asinh(estop) / N;
W3  = Wstop / cosh(acosh(estop)/N);
f3a = (fs/pi) * atan((sqrt(W3^2 - c^2 + 1) - W3)/(c+1));
f3b = (fs/pi) * atan((sqrt(W3^2 - c^2 + 1) + W3)/(c+1));
P = [Wpass,Wstop,Wsa,Wsb,c,epass,estop,Nex,N,f3a,f3b,a];
W0 = sinh(a) / Wstop;                         % reciprocal of text

if r==1,
    G = 1 / (1 + W0);
    a1 = -2 * c * W0 / (1 + W0);
    a2 = -(1 - W0) / (1 + W0);
    A(1,:) = [1, a1, a2, 0, 0];
    B(1,:) = G * [1, 0, -1, 0, 0];
else
    A(1,:) = [1, 0, 0, 0, 0];
    B(1,:) = [1, 0, 0, 0, 0];
end

for i=1:K,
    th = pi * (N - 1 + 2 * i) / (2 * N);
    Wi = sin(th) / Wstop;                       % reciprocal of text
    D = 1 - 2 * W0 * cos(th) + W0^2 + Wi^2;
    G = (1 + Wi^2) / D;
    b1 = - 4 * c * Wi^2 / (1 + Wi^2);
    b2 = 2 * (Wi^2 * (2*c*c+1) - 1) / (1 + Wi^2);
    a1 = 4 * c * (W0 * cos(th) - W0^2 - Wi^2) / D;
    a2 = 2 * ((2*c*c + 1)*(W0^2 + Wi^2) - 1) / D;
    a3 = - 4 * c * (W0 * cos(th) + W0^2 + Wi^2) / D;
    a4 = (1 + 2 * W0 * cos(th) + W0^2 + Wi^2) / D;
    A(i+1,:) = [1, a1, a2, a3, a4];
    B(i+1,:) = G * [1, b1, b2, b1, 1];
end
