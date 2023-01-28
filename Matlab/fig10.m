% fig10.m - Butterworth example
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

type = 0; 

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

G0 = 0; Gb=3;
G1 = 9;  GB1 = G1-Gb;   % GB is 3-dB below peak gain
G2 = 12; GB2 = G2-Gb; 
G3 = -6; GB3 = G3+Gb;   % GB is 3-dB above cut gain
G4 = 6;  GB4 = G4-Gb; 

[B1,A1,Bh1,Ah1] = hpeq(N, G0, G1, GB1, w1, Dw1, type); 
[B2,A2,Bh2,Ah2] = hpeq(N, G0, G2, GB2, w2, Dw2, type);
[B3,A3,Bh3,Ah3] = hpeq(N, G0, G3, GB3, w3, Dw3, type);
[B4,A4,Bh4,Ah4] = hpeq(N, G0, G4, GB4, w4, Dw4, type);  

% biquad version:

[Bq1,Aq1] = hpeq(1, G0, G1, GB1, w1, Dw1, 0);    % individual biquads
[Bq2,Aq2] = hpeq(1, G0, G2, GB2, w2, Dw2, 0);
[Bq3,Aq3] = hpeq(1, G0, G3, GB3, w3, Dw3, 0);
[Bq4,Aq4] = hpeq(1, G0, G4, GB4, w4, Dw4, 0);  


f = linspace(0,20,501); w = 2*pi*f/fs;

% cascaded equalizers - each shifted to its own center frequency
H = fresp(Bh1,Ah1,w,w1) .* fresp(Bh2,Ah2,w,w2) .* fresp(Bh3,Ah3,w,w3) .* fresp(Bh4,Ah4,w,w4);

% equivalently, use the cascaded 4-th order sections:
% H = fresp([B1; B2; B3; B4],[A1; A2; A3; A4], w);

% cascaded biquads 
Hbquad = fresp([Bq1; Bq2; Bq3; Bq4],[Aq1; Aq2; Aq3; Aq4], w);

H = 20*log10(abs(H));
Hbquad = 20*log10(abs(Hbquad));

% individual biquads
Hb1 = 20*log10(abs(fresp(Bq1,Aq1,w)));
Hb2 = 20*log10(abs(fresp(Bq2,Aq2,w)));
Hb3 = 20*log10(abs(fresp(Bq3,Aq3,w)));
Hb4 = 20*log10(abs(fresp(Bq4,Aq4,w)));

figure;

plot(f,H,'b-', f,Hbquad,'--r', f,Hb1,'m:', f,Hb2,'g:', f,Hb3,'c:', f,Hb4,'k:');
hold on;

legend('cascaded Butterworth','cascaded biquads', 'individual biquads',1);

ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
grid;

plot([f12, f21,f2,f22, f31,f3,f32, f41], [GB1, GB2,G2,GB2, GB3,G3,GB3, GB4], '.', 'markersize',18);
hold off;

xlabel('{\it f}  (kHz)');
ylabel('dB');
title('{\it N}=4, Butterworth');

print -depsc fig10a.eps


% ------------------------------------------------------------------------

N = 5;

[B1,A1,Bh1,Ah1] = hpeq(N, G0, G1, GB1, w1, Dw1, type); 
[B2,A2,Bh2,Ah2] = hpeq(N, G0, G2, GB2, w2, Dw2, type);
[B3,A3,Bh3,Ah3] = hpeq(N, G0, G3, GB3, w3, Dw3, type);
[B4,A4,Bh4,Ah4] = hpeq(N, G0, G4, GB4, w4, Dw4, type);  

% biquad version:

[Bq1,Aq1] = hpeq(1, G0, G1, GB1, w1, Dw1, 0);  
[Bq2,Aq2] = hpeq(1, G0, G2, GB2, w2, Dw2, 0);
[Bq3,Aq3] = hpeq(1, G0, G3, GB3, w3, Dw3, 0);
[Bq4,Aq4] = hpeq(1, G0, G4, GB4, w4, Dw4, 0);  

% cascaded equalizers - each shifted to its own center frequency
H = fresp(Bh1,Ah1,w,w1) .* fresp(Bh2,Ah2,w,w2) .* fresp(Bh3,Ah3,w,w3) .* fresp(Bh4,Ah4,w,w4);

% alternatively, cascade the fourth-order sections, which are already frequency-shifted
% H = fresp([B1; B2; B3; B4],[A1; A2; A3; A4], w);

% cascaded biquads - biquads are already frequency-shifted
Hbquad = fresp([Bq1; Bq2; Bq3; Bq4],[Aq1; Aq2; Aq3; Aq4], w);

H = 20*log10(abs(H));
Hbquad = 20*log10(abs(Hbquad));

figure;

plot(f,H,'b-', f,Hbquad,'--r');
hold on;

ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
grid;

plot([f12, f21,f2,f22, f31,f3,f32, f41], [GB1, GB2,G2,GB2, GB3,G3,GB3, GB4], '.', 'markersize',18);
hold off;

xlabel('{\it f}  (kHz)');
ylabel('dB');
title('{\it N}=5, Butterworth');

print -depsc fig10b.eps

