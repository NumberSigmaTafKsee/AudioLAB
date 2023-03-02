% fig14.m - redesigned with common 3-dB widths
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

fs=40; 
f1 = 0;  Df1 = 1;  w1 = 2*pi*f1/fs; Dwb1 = 2*pi*Df1/fs;          
f2 = 4;  Df2 = 2;  w2 = 2*pi*f2/fs; Dwb2 = 2*pi*Df2/fs;
f3 = 9;  Df3 = 2;  w3 = 2*pi*f3/fs; Dwb3 = 2*pi*Df3/fs;
f4 = 20; Df4 = 4;  w4 = 2*pi*f4/fs; Dwb4 = 2*pi*Df4/fs;

[w11,w12] = bandedge(w1,Dwb1); fb11 = w11*fs/2/pi; fb12 = w12*fs/2/pi;   % 3-dB bandedges
[w21,w22] = bandedge(w2,Dwb2); fb21 = w21*fs/2/pi; fb22 = w22*fs/2/pi; 
[w31,w32] = bandedge(w3,Dwb3); fb31 = w31*fs/2/pi; fb32 = w32*fs/2/pi; 
[w41,w42] = bandedge(w4,Dwb4); fb41 = w41*fs/2/pi; fb42 = w42*fs/2/pi; 

f = linspace(0,20,501); w = 2*pi*f/fs;

N = 4;

% ---------------------------- Butterworth --------------------------

type = 0;  

G0 = 0; Gb=3;

G1 = 9;  GB1 = 6;   % 3-dB below peak gain
G2 = 12; GB2 = 9; 
G3 = -6; GB3 = -3; 
G4 = 6;  GB4 = 3; 

[B1,A1] = hpeq(N, G0, G1, GB1, w1, Dwb1, type); 
[B2,A2] = hpeq(N, G0, G2, GB2, w2, Dwb2, type);
[B3,A3] = hpeq(N, G0, G3, GB3, w3, Dwb3, type);
[B4,A4] = hpeq(N, G0, G4, GB4, w4, Dwb4, type);  

B = [B1; B2; B3; B4];
A = [A1; A2; A3; A4];

H = 20*log10(abs(fresp(B,A,w)));

figure; 

plot(f,H,'b-');
hold on;

ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
grid;

plot([fb12, fb21,f2,fb22, fb31,f3,fb32, fb41], [GB1, GB2,G2,GB2, GB3,G3,GB3, GB4], 'k.', 'markersize',18);
hold off;

xlabel('{\it f}  (kHz)');
ylabel('dB');
title('{\it N} = 4, Butterworth');

print -deps fig14a.eps

% ---------------------------- Chebyshev-1 --------------------------

type = 1;  

G0 = 0; Gb=0.01;  gb = 3;

G1 = 9;  GB1 = 8.99;  Gb1 = 6;
G2 = 12; GB2 = 11.99; Gb2 = 9;
G3 = -6; GB3 = -5.99; Gb3 = -3;
G4 = 6;  GB4 = 5.99;  Gb4 = 3;

Dw1 = hpeqbw(N, G0, G1, GB1, Gb1, Dwb1, type);    % Dw1 is bandwidth at level GB1
Dw2 = hpeqbw(N, G0, G2, GB2, Gb2, Dwb2, type);
Dw3 = hpeqbw(N, G0, G3, GB3, Gb3, Dwb3, type);
Dw4 = hpeqbw(N, G0, G4, GB4, Gb4, Dwb4, type);

[w11,w12] = bandedge(w1,Dw1); f11 = w11*fs/2/pi; f12 = w12*fs/2/pi;   % bandedges at level GB1
[w21,w22] = bandedge(w2,Dw2); f21 = w21*fs/2/pi; f22 = w22*fs/2/pi; 
[w31,w32] = bandedge(w3,Dw3); f31 = w31*fs/2/pi; f32 = w32*fs/2/pi; 
[w41,w42] = bandedge(w4,Dw4); f41 = w41*fs/2/pi; f42 = w42*fs/2/pi; 

[B1,A1] = hpeq(N, G0, G1, GB1, w1, Dw1, type); 
[B2,A2] = hpeq(N, G0, G2, GB2, w2, Dw2, type);
[B3,A3] = hpeq(N, G0, G3, GB3, w3, Dw3, type);
[B4,A4] = hpeq(N, G0, G4, GB4, w4, Dw4, type);  

B = [B1; B2; B3; B4];
A = [A1; A2; A3; A4];

H = 20*log10(abs(fresp(B,A,w)));

figure; 

plot(f,H,'b-');
hold on;

ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
grid;

plot([fb12, fb21,f2,fb22, fb31,f3,fb32, fb41], [Gb1, Gb2,G2,Gb2, Gb3,G3,Gb3, Gb4], 'k.', 'markersize',18);
hold off;

xlabel('{\it f}  (kHz)');
ylabel('dB');
title('{\it N} = 4, Chebyshev-1');

print -deps fig14b.eps

% ---------------------------- Chebyshev-2 --------------------------

type = 2;  

G0 = 0; Gb=0.01;  gb = 3;

G1 = 9;  GB1 = 0.01;  Gb1 = 6;
G2 = 12; GB2 = 0.01;  Gb2 = 9;
G3 = -6; GB3 = -0.01; Gb3 = -3;
G4 = 6;  GB4 = 0.01;  Gb4 = 3;

Dw1 = hpeqbw(N, G0, G1, GB1, Gb1, Dwb1, type);
Dw2 = hpeqbw(N, G0, G2, GB2, Gb2, Dwb2, type);
Dw3 = hpeqbw(N, G0, G3, GB3, Gb3, Dwb3, type);
Dw4 = hpeqbw(N, G0, G4, GB4, Gb4, Dwb4, type);

[w11,w12] = bandedge(w1,Dw1); f11 = w11*fs/2/pi; f12 = w12*fs/2/pi; 
[w21,w22] = bandedge(w2,Dw2); f21 = w21*fs/2/pi; f22 = w22*fs/2/pi; 
[w31,w32] = bandedge(w3,Dw3); f31 = w31*fs/2/pi; f32 = w32*fs/2/pi; 
[w41,w42] = bandedge(w4,Dw4); f41 = w41*fs/2/pi; f42 = w42*fs/2/pi; 

[B1,A1] = hpeq(N, G0, G1, GB1, w1, Dw1, type); 
[B2,A2] = hpeq(N, G0, G2, GB2, w2, Dw2, type);
[B3,A3] = hpeq(N, G0, G3, GB3, w3, Dw3, type);
[B4,A4] = hpeq(N, G0, G4, GB4, w4, Dw4, type);  

B = [B1; B2; B3; B4];
A = [A1; A2; A3; A4];

H = 20*log10(abs(fresp(B,A,w)));

figure; 

plot(f,H,'b-');
hold on;

ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
grid;

plot([fb12, fb21,f2,fb22, fb31,f3,fb32, fb41], [Gb1, Gb2,G2,Gb2, Gb3,G3,Gb3, Gb4], 'k.', 'markersize',18);
hold off;

xlabel('{\it f}  (kHz)');
ylabel('dB');
title('{\it N} = 4, Chebyshev-2');

print -deps fig14c.eps


% ---------------------------- Elliptic --------------------------

type = 3;  

G0 = 0; Gb=0.01;  gb = 3;

G1 = 9;  GB1 = 8.99;  Gs1 = 0.01;  Gb1 = 6;
G2 = 12; GB2 = 11.99; Gs2 = 0.01;  Gb2 = 9;
G3 = -6; GB3 = -5.99; Gs3 = -0.01; Gb3 = -3;
G4 = 6;  GB4 = 5.99;  Gs4 = 0.01;  Gb4 = 3;

Dw1 = hpeqbw(N, G0, G1, GB1, Gb1, Dwb1, type, Gs1);
Dw2 = hpeqbw(N, G0, G2, GB2, Gb2, Dwb2, type, Gs2);
Dw3 = hpeqbw(N, G0, G3, GB3, Gb3, Dwb3, type, Gs3);
Dw4 = hpeqbw(N, G0, G4, GB4, Gb4, Dwb4, type, Gs4);

[w11,w12] = bandedge(w1,Dw1); f11 = w11*fs/2/pi; f12 = w12*fs/2/pi; 
[w21,w22] = bandedge(w2,Dw2); f21 = w21*fs/2/pi; f22 = w22*fs/2/pi; 
[w31,w32] = bandedge(w3,Dw3); f31 = w31*fs/2/pi; f32 = w32*fs/2/pi; 
[w41,w42] = bandedge(w4,Dw4); f41 = w41*fs/2/pi; f42 = w42*fs/2/pi; 

[B1,A1,B1h,A1h] = hpeq(N, G0, G1, GB1, w1, Dw1, type, Gs1);  
[B2,A2,B2h,A2h] = hpeq(N, G0, G2, GB2, w2, Dw2, type, Gs2); 
[B3,A3,B3h,A3h] = hpeq(N, G0, G3, GB3, w3, Dw3, type, Gs3); 
[B4,A4,B4h,A4h] = hpeq(N, G0, G4, GB4, w4, Dw4, type, Gs4); 

B = [B1; B2; B3; B4];
A = [A1; A2; A3; A4];

H = 20*log10(abs(fresp(B,A,w)));

figure; 

plot(f,H,'b-');
hold on;

ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
grid;

plot([fb12, fb21,f2,fb22, fb31,f3,fb32, fb41], [Gb1, Gb2,G2,Gb2, Gb3,G3,Gb3, Gb4], 'k.', 'markersize',18);
hold off;

xlabel('{\it f}  (kHz)');
ylabel('dB');
title('{\it N} = 4, Elliptic');

print -deps fig14d.eps


