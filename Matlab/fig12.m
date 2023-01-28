% fig12.m - Chebyshev-2 example
%
% -------------------------------------------------------------------------
% Copyright (c) 2005 by Sophocles J. Orfanidis
% 
% Address: Sophocles J. Orfanidis                       
%          ECE Department, Rutgers University          
%          94 Brett Road, Piscataway, NJ 08854-8058, USA
%
% Email:   orfanidi@ece.rutgers.edu
% Date:    June 15, 2005
% 
% Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
%
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% 
% tested with MATLAB R11.1 and R14
% -------------------------------------------------------------------------

clear all;

N = 4;
fs=40; 
f1 = 0;  Df1 = 1;  w1 = 2*pi*f1/fs; Dw1 = 2*pi*Df1/fs;
f2 = 4;  Df2 = 2;  w2 = 2*pi*f2/fs; Dw2 = 2*pi*Df2/fs;
f3 = 9;  Df3 = 2;  w3 = 2*pi*f3/fs; Dw3 = 2*pi*Df3/fs;
f4 = 20; Df4 = 4;  w4 = 2*pi*f4/fs; Dw4 = 2*pi*Df4/fs;

[w11,w12] = bandedge(w1,Dw1); f11 = w11*fs/2/pi; f12 = w12*fs/2/pi; 
[w21,w22] = bandedge(w2,Dw2); f21 = w21*fs/2/pi; f22 = w22*fs/2/pi; 
[w31,w32] = bandedge(w3,Dw3); f31 = w31*fs/2/pi; f32 = w32*fs/2/pi; 
[w41,w42] = bandedge(w4,Dw4); f41 = w41*fs/2/pi; f42 = w42*fs/2/pi; 

G0 = 0;  Gb=0.01;
G1 = 9;  GB1 = G0+Gb; 
G2 = 12; GB2 = G0+Gb; 
G3 = -6; GB3 = G0-Gb; 
G4 = 6;  GB4 = G0+Gb; 

type = 2; 

[B1,A1] = hpeq(N, G0, G1, GB1, w1, Dw1, type); 
[B2,A2] = hpeq(N, G0, G2, GB2, w2, Dw2, type);
[B3,A3] = hpeq(N, G0, G3, GB3, w3, Dw3, type);
[B4,A4] = hpeq(N, G0, G4, GB4, w4, Dw4, type);  

B = [B1; B2; B3; B4];
A = [A1; A2; A3; A4];

f = linspace(0,20,501); w = 2*pi*f/fs;

H = 20*log10(abs(fresp(B,A,w)));

figure; 

plot(f,H,'b-');
hold on;

ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
grid;

plot([f12, f21,f2,f22, f31,f3,f32, f41], [GB1, GB2,G2,GB2, GB3,G3,GB3, GB4], '.', 'markersize',18);
hold off;

xlabel('{\it f}  (kHz)');
ylabel('dB');
title('{\it N}=4, Chebyshev, type-2');
print -depsc fig12a.eps


N = 5;

[B1,A1] = hpeq(N, G0, G1, GB1, w1, Dw1, type); 
[B2,A2] = hpeq(N, G0, G2, GB2, w2, Dw2, type);
[B3,A3] = hpeq(N, G0, G3, GB3, w3, Dw3, type);
[B4,A4] = hpeq(N, G0, G4, GB4, w4, Dw4, type);  

B = [B1; B2; B3; B4];
A = [A1; A2; A3; A4];

H = 20*log10(abs(fresp(B,A,w)));

figure; 

%plot(f,H,'b-', f,Hbquad,'--r');
plot(f,H,'b-');
hold on;

ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
grid;

plot([f12, f21,f2,f22, f31,f3,f32, f41], [GB1, GB2,G2,GB2, GB3,G3,GB3, GB4], '.', 'markersize',18);
hold off;

xlabel('{\it f}  (kHz)');
ylabel('dB');
title('{\it N}=5, Chebyshev, type-2');
print -depsc fig12b.eps



