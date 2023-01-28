% hpeqex2.m - example output from hpeq.m
%
% high-order designs with given 3-dB width

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

fs = 40;
tol = eps;

disp(' ');
disp('Equalizer designs based on 3-dB width')
disp('-------------------------------------');
disp('sampling rate is 40 kHz');
disp('all gains must be entered in dB (enter -Inf to get 0 in absolute units)');
disp('for plotting purposes, limit peak gain G to [-15,15] dB, unless G=-Inf');

disp(' ');
disp('enter the design parameters as follows, with f''s in kHz (include the brackets):');
disp(' ');
disp('   [N, G0, G, GB, Gb, f0, Dfb, type]                 (type=0,1,2)      (Dfb is bandwidth at level Gb)');
disp('   [N, G0, G, GB, Gb, f0, Dfb, type, Gs]             (type=3)');
disp('   [m, G0, G, GB, Gb, f1, f2,  f1b,  f2b, type]      (N to be determined, m=1,2, type=0,1,2)');
disp('   [m, G0, G, GB, Gb, f1, f2,  f1b,  f2b, type, Gs]  (N to be determined, m=1,2, type=3)');
disp('                                                     (m=1 matches Dw exactly, m=2 matches Dwb)');

%disp(' ');
disp('EQ and shelving examples with fixed order:');
disp('------------------------------------------');
disp('[5, 0, 12, 9,    9, 6, 4, 0],      [5, 0, 12, 9,    9, 0, 4, 0],      [5, 0, 12, 9,    20, 4, 0]         = Butter  EQ,LP,HP');
disp('[5, 0, 12, 11.5, 9, 6, 4, 1],      [5, 0, 12, 11.5, 9, 0, 4, 1],      [5, 0, 12, 11.5, 20, 4, 1]         = Cheby-1 EQ,LP,HP');
disp('[5, 0, 12, 0.5,  9, 6, 4, 2],      [5, 0, 12, 0.5,  9, 0, 4, 2],      [5, 0, 12, 0.5,  20, 4, 2]         = Cheby-2 EQ,LP,HP');
disp('[5, 0, 12, 11.5, 9, 6, 4, 3, 0.5], [5, 0, 12, 11.5, 9, 0, 4, 3, 0.5], [5, 0, 12, 11.5, 9, 20, 4, 3, 0.5] = Ellip   EQ,LP,HP');

disp(' ');
disp('EQ examples with unknown order (set f1=f1b=0 for LP, or, f2=f2b=20 for HP):');
disp('---------------------------------------------------------------------------');
disp('[1, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 0],      [2, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 0]      = Butter,  m=1,2');
disp('[1, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 1],      [2, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 1]      = Cheby-1, m=1,2');
disp('[1, 0, 12, 0.1,  9, 5.5, 8.5, 6,   8,   2],      [2, 0, 12, 0.1,  9, 5.5, 8.5, 6,   8,   2]      = Cheby-2, m=1,2');
disp('[1, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 3, 0.1], [2, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 3, 0.1] = Ellip,   m=1,2');

disp(' ');
disp('Ordinary BP, BS with fixed order (set f0=0,20 for LP,HP):');
disp('--------------------------------------------------------- ');
disp('[5, -inf, 0, -1,  -3, 6, 4, 0],      [5, 0, -inf, -20, -3, 6, 4, 0]     = Butter  EQ,LP,HP');
disp('[5, -inf, 0, -1,  -3, 6, 4, 1],      [5, 0, -inf, -20, -3, 6, 4, 1]     = Cheby-1 EQ,LP,HP');
disp('[5, -inf, 0, -20, -3, 6, 4, 2],      [5, 0, -inf, -1,  -3, 6, 4, 2]     = Cheby-2 EQ,LP,HP');
disp('[5, -inf, 0, -1,  -3, 6, 4, 3, -20], [5, 0, -inf, -20, -3, 6, 4, 3, -1] = Ellip   EQ,LP,HP');

disp(' ');
disp('Ordinary BP, BS with unknown order (set f1=f1b=0 for LP, or, f2=f2b=20 for HP):');
disp('------------------------------------------------------------------------------- ');
disp('[1, -inf, 0, -1,  -3, 6,   8,   5.5, 8.5, 0],      [2, -inf, 0, -1,  -3, 6,   8,   5.5, 8.5, 0]      = Butter,  m=1,2');
disp('[1, -inf, 0, -1,  -3, 6,   8,   5.5, 8.5, 1],      [2, -inf, 0, -1,  -3, 6,   8,   5.5, 8.5, 1]      = Cheby-1, m=1,2');
disp('[1, -inf, 0, -20, -3, 5.5, 8.5, 6,   8,   2],      [2, -inf, 0, -20, -3, 5.5, 8.5, 6,   8,   2]      = Cheby-2, m=1,2');
disp('[1, -inf, 0, -1,  -3, 6,   8,   5.5, 8.5, 3, -20], [2, -inf, 0, -1,  -3, 6,   8,   5.5, 8.5, 3, -20] = Ellip,   m=1,2');
disp(' ');

P = input('enter design parameters (or, copy one of the above, 0 to exit) = ');

if P==0; return; end

switch length(P)
   case 8,
      N = P(1); G0 = P(2); G = P(3); GB = P(4); Gb = P(5); f0 = P(6); Dfb = P(7); type = P(8); 
      w0 = 2*pi*f0/fs; Dwb = 2*pi*Dfb/fs;
   case 9,
      N = P(1); G0 = P(2); G = P(3); GB = P(4); Gb = P(5); f0 = P(6); Dfb = P(7); type = P(8); Gs = P(9);
      w0 = 2*pi*f0/fs; Dwb = 2*pi*Dfb/fs;
   case 10,
      m = P(1); G0 = P(2); G = P(3); GB = P(4); Gb = P(5); 
      f1 = P(6); f2 = P(7); f1b = P(8); f2b = P(9); type = P(10);
      w1 = 2*pi*f1/fs; w2 = 2*pi*f2/fs; w1b = 2*pi*f1b/fs; w2b = 2*pi*f2b/fs;
      [w0,Dw] = bandedge(w1,w2,1);  [w0b,Dwb] = bandedge(w1b,w2b,1);
      f0 = w0 * fs/2/pi; f0b = w0b * fs/2/pi;
   case 11,
      m = P(1); G0 = P(2); G = P(3); GB = P(4); Gb = P(5); 
      f1 = P(6); f2 = P(7); f1b = P(8); f2b = P(9); type = P(10); Gs = P(11);
      w1 = 2*pi*f1/fs; w2 = 2*pi*f2/fs; w1b = 2*pi*f1b/fs; w2b = 2*pi*f2b/fs;
      [w0,Dw] = bandedge(w1,w2,1);  [w0b,Dwb] = bandedge(w1b,w2b,1);
      f0 = w0 * fs/2/pi; f0b = w0b * fs/2/pi;
end

switch length(P)
   case 8,
      Dw = hpeqbw(N,G0,G,GB,Gb,Dwb,type);
      [w1,w2] = bandedge(w0,Dw);    f1  = w1 * fs/2/pi;  f2  = w2 * fs/2/pi;
      [w1b,w2b] = bandedge(w0,Dwb); f1b = w1b * fs/2/pi; f2b = w2b * fs/2/pi;
   case 9,
      Dw = hpeqbw(N,G0,G,GB,Gb,Dwb,3,Gs);
      [w1,w2] = bandedge(w0,Dw);    f1  = w1 * fs/2/pi;  f2  = w2 * fs/2/pi;
      [w1b,w2b] = bandedge(w0,Dwb); f1b = w1b * fs/2/pi; f2b = w2b * fs/2/pi;
   case 10,
      N = ceil(hpeqord(G0,G,GB,Gb,Dw,Dwb,type));
      if m==2, 
         Dw = hpeqbw(N,G0,G,GB,Gb,Dwb,type); w0 = w0b; f0 = w0 * fs/2/pi;
      end
   case 11,
      N = ceil(hpeqord(G0,G,GB,Gs,Dw,Dwb,3,Gb));
      if m==2, 
         Dw = hpeqbw(N,G0,G,GB,Gb,Dwb,3,Gs); w0 = w0b; f0 = w0 * fs/2/pi;
      end
end

r = mod(N,2);   

g0 = 10^(G0/20); g = 10^(G/20); gB = 10^(GB/20); gb = 10^(Gb/20);    
if length(P)==9 | length(P)==11, gs = 10^(Gs/20); end

switch type,                                        % H0 and Hinf gains in absolute units
   case 0,
      h0 = g; hinf = g0;
   case 1,
      h0 = g^r * gB^(1-r); hinf = g0;
   case 2,
      h0 = g; hinf = g0^r * gB^(1-r);
   case 3,
      h0 = g^r * gB^(1-r); hinf = g0^r * gs^(1-r);
end

H0 = 20*log10(h0); Hinf = 20*log10(hinf);           % H0 and Hinf gains in dB     

if type==3,
   [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,3,Gs,tol)    
else 
   [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,type)
end

figure;

f = linspace(0,20,1001); w = 2*pi*f/fs;

H = abs(fresp(B,A,w));

if (g0~=0) & (g~=0),             % equalizer + shelving filters

  H = 20*log10(H);
  plot(f,H,'b-', [0,f1b,f1,f0,f2,f2b,20], [Hinf,Gb,GB,H0,GB,Gb,Hinf],'r.');  
  ylim([-15,15]); ytick(-15:3:15); ylabel('dB');
  xtick(0:2:20); xlabel('{\it f}  (kHz)'); grid;

else                              % ordinary bandpass, bandstop, lowpass, highpass 

  plot(f,H,'b-', [0,f1b,f1,f0,f2,f2b,20], [hinf,gb,gB,h0,gB,gb,hinf], 'r.');  
  ylim([0,1.05]); ytick(0:0.1:1); ylabel('absolute units');
  xtick(0:2:20); xlabel('{\it f}  (kHz)'); grid;

end

