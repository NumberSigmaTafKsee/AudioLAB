% lhbutt.m - lowpass/highpass Butterworth digital filter design
%
% [A, B, P] = lhbutt(s, fs, fpass, fstop, Apass, Astop)
%
% s = 1, -1 = lowpass, highpass
% design parameters:
% P = [Wpass, Wstop, epass, estop, Nex, N, Astop, W0, f0];
% A, B are Kx3 matrices of cascade second-order sections

function [A, B, P] = lhbutt(s, fs, fpass, fstop, Apass, Astop)

Wpass = tan(pi * fpass / fs);  Wpass = Wpass^s;         % cot() for HP
Wstop = tan(pi * fstop / fs);  Wstop = Wstop^s;

epass = sqrt(10^(Apass/10) - 1); 
estop = sqrt(10^(Astop/10) - 1);

Nex = log(estop/epass) / log(Wstop/Wpass);
N = ceil(Nex);  r = rem(N,2);  K = (N - r) / 2;         % K = no. sections

W0 = Wpass * (epass^(-1/N));
Astop = 10 * log10(1 + (Wstop/W0)^(2*N));               % actual Astop
f0 = (fs/pi) * atan(W0^s);                              % 3-dB freq. in Hz
P = [Wpass, Wstop, epass, estop, Nex, N, Astop, W0, f0];

if r==1,                                                % N = odd
    G = W0 / (1 + W0);                                  % 1st order section
    B(1,:) = G * [1, s, 0];
    A(1,:) = [1, s*(2*G-1), 0];
else                                                    % N = even
    B(1,:) = [1, 0, 0];
    A(1,:) = [1, 0, 0];
end

for i=1:K,
    th = pi * (N - 1 + 2 * i) / (2 * N);
    D = 1 - 2 * W0 * cos(th) + W0^2;
    G = W0^2 / D;
    a1 = 2 * (W0^2 - 1) / D;
    a2 = (1 + 2 * W0 * cos(th) + W0^2) / D;
    B(i+1,:) = G * [1, 2*s, 1];
    A(i+1,:) = [1, s*a1, a2];
end
