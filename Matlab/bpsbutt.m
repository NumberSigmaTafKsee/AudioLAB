% bpsbutt.m - bandpass/bandstop Butterworth digital filter design
%
% [A, B, P] = bpsbutt(s, fs, fpa, fpb, fsa, fsb, Apass, Astop)
%
% s = 1, -1 = bandpass, bandstop
% design parameters:
% P = [Wpass,Wstop,Wsa,Wsb,c,fc,epass,estop,Nex,N,Astop,W0,f0a,f0b];
% A, B are Kx5 matrices of cascade of fourth-order sections

function [A, B, P] = bpsbutt(s, fs, fpa, fpb, fsa, fsb, Apass, Astop)

c = sin(2*pi*(fpa + fpb)/fs) / (sin(2*pi*fpa/fs) + sin(2*pi*fpb/fs));
fc = 0.5 * (fs/pi) * acos(c);
Wpass = (abs((c - cos(2*pi*fpb/fs)) / sin(2*pi*fpb/fs)))^s;
Wsa = (c - cos(2*pi*fsa/fs)) / sin(2*pi*fsa/fs);  Wsa = Wsa^s;
Wsb = (c - cos(2*pi*fsb/fs)) / sin(2*pi*fsb/fs);  Wsb = Wsb^s;
Wstop = min(abs(Wsa), abs(Wsb));

epass = sqrt(10^(Apass/10) - 1);
estop = sqrt(10^(Astop/10) - 1);

Nex = log(estop/epass) / log(Wstop/Wpass);
N = ceil(Nex);  r = rem(N,2);  K = (N - r) / 2;

W0 = Wpass * (epass^(-1/N));  W0s = W0^s;
Astop = 10 * log10(1 + (Wstop/W0)^(2*N));
f0a = (fs/pi) * atan((sqrt(W0s^2 - c^2 + 1) - W0s)/(c+1));
f0b = (fs/pi) * atan((sqrt(W0s^2 - c^2 + 1) + W0s)/(c+1));
P = [Wpass,Wstop,Wsa,Wsb,c,fc,epass,estop,Nex,N,Astop,W0,f0a,f0b];

if r==1,
    G = W0 / (1 + W0);
    a1 = -2 * c / (1 + W0s);
    a2 = (1 - W0s) / (1 + W0s);
    A(1,:) = [1, a1, a2, 0, 0];
    B(1,:) = G * [1, (s-1)*c, -s, 0, 0];
else
    A(1,:) = [1, 0, 0, 0, 0];
    B(1,:) = [1, 0, 0, 0, 0];
end

for i=1:K,
    th = pi * (N - 1 + 2 * i) / (2 * N);
    D = 1 - 2 * W0s * cos(th) + W0s^2;
    G = W0^2 / (1 - 2 * W0 * cos(th) + W0^2);
    a1 = 4 * c * (W0s * cos(th) - 1) / D;
    a2 = 2 * (2*c^2 + 1 - W0s^2) / D;
    a3 = - 4 * c * (W0s * cos(th) + 1) / D;
    a4 = (1 + 2 * W0s * cos(th) + W0s^2) / D;
    A(i+1,:) = [1, a1, a2, a3, a4];
    B(i+1,:) = G * conv([1, (s-1)*c, -s], [1, (s-1)*c, -s]);
end
