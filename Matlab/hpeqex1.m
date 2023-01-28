% hpeqex1.m - example output from hpeq.m
%
% EQ and LP,HP shelving examples with fixed or unknown order
% Ordinary LP,HP,BP,BS examples with fixed or unknown order

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
disp('sampling rate is 40 kHz');
disp('all gains must be entered in dB (enter -Inf to get 0 in absolute units)');
disp('for plotting purposes, limit peak gain G to [-15,15] dB, unless G=-Inf');
disp(' ');
disp('enter the design parameters as follows, with f''s in kHz (include the brackets):');
disp(' ');
disp('   [N,  G0, G,  GB, f0, Df, type]            (for type=0,1,2) ');
disp('   [N,  G0, G,  GB, f0, Df, type, Gs]        (for type=3)');
disp('   [G0, G,  GB, Gs, f1, f2, f1s, f2s, type]  (with N to be determined)');

disp(' ');
disp('EQ and shelving examples with fixed or unknown order:');
disp('----------------------------------------------------- ');
disp('[5, 0, 12,    9, 6, 4, 0],    [5, 0, 12,    9, 0, 4, 0],    [5, 0, 12,    9, 20, 4, 0]    = Butter  EQ,LP,HP');
disp('[5, 0, 12, 11.5, 6, 4, 1],    [5, 0, 12, 11.5, 0, 4, 1],    [5, 0, 12, 11.5, 20, 4, 1]    = Cheby-1 EQ,LP,HP');
disp('[5, 0, 12,  0.5, 6, 4, 2],    [5, 0, 12,  0.5, 0, 4, 2],    [5, 0, 12,  0.5, 20, 4, 2]    = Cheby-2 EQ,LP,HP');
disp('[5, 0, 12, 11.5, 6, 4, 3, 1], [5, 0, 12, 11.5, 0, 4, 3, 1], [5, 0, 12, 11.5, 20, 4, 3, 1] = Ellip   EQ,LP,HP');
disp(' ');
disp('[0, 12, 11.5,  0.5, 6, 8, 5, 9, 0], [0, 12, 11.5,  0.5, 0, 4, 0, 5, 0], [0, 12, 11.5,  0.5, 16, 20, 15, 20, 0] = Butter');
disp('[0, 12, 11.5,  0.5, 6, 8, 5, 9, 1], [0, 12, 11.5,  0.5, 0, 4, 0, 5, 1], [0, 12, 11.5,  0.5, 16, 20, 15, 20, 1] = Cheby-1');
disp('[0, 12,  0.5, 11.5, 5, 9, 6, 8, 2], [0, 12,  0.5, 11.5, 0, 5, 0, 4, 2], [0, 12,  0.5, 11.5, 15, 20, 16, 20, 2] = Cheby-2');
disp('[0, 12, 11.5,  0.5, 6, 8, 5, 9, 3], [0, 12, 11.5,  0.5, 0, 4, 0, 5, 3], [0, 12, 11.5,  0.5, 16, 20, 15, 20, 3] = Ellip');

disp(' ');
disp('Ordinary BP, BS with fixed or unknown order:');
disp('-------------------------------------------- ');
disp('[5, -inf, 0,  -3, 6, 4, 0],       [5, 0, -inf,  -3, 6, 4, 0]     = Butter');
disp('[5, -inf, 0,  -1, 6, 4, 1],       [5, 0, -inf, -20, 6, 4, 1]     = Cheby-1');
disp('[5, -inf, 0, -20, 6, 4, 2],       [5, 0, -inf,  -1, 6, 4, 2]     = Cheby-2');
disp('[5, -inf, 0,  -1, 6, 4, 3, -20],  [5, 0, -inf, -20, 6, 4, 3, -1] = Ellip');
disp(' ');
disp('[-inf, 0,  -1, -20, 6, 8, 5, 9, t],  [0, -inf, -20,  -1, 6, 8, 5, 9, t] = Butter, Cheby-1, Ellip, t=0,1,3');
disp('[-inf, 0, -20,  -1, 5, 9, 6, 8, t],  [0, -inf,  -1, -20, 5, 9, 6, 8, t] = Cheby-2, t=2');

disp(' ');
disp('Ordinary LP, HP with fixed or unknown order:');
disp('-------------------------------------------- ');
disp('[5, -inf, 0,  -3, 0, 4, 0],       [5, -inf, 0,  -3, 20, 4, 0]      = Butter');
disp('[5, -inf, 0,  -1, 0, 4, 1],       [5, -inf, 0,  -1, 20, 4, 1]      = Cheby-1');
disp('[5, -inf, 0, -20, 0, 4, 2],       [5, -inf, 0, -20, 20, 4, 2]      = Cheby-2');
disp('[5, -inf, 0,  -1, 0, 4, 3, -20],  [5, -inf, 0,  -1, 20, 4, 3, -20] = Ellip');
disp(' ');
disp('[-inf, 0,  -1, -20, 0, 4, 0, 5, t],  [-inf, 0,  -1, -20, 16, 20, 15, 20, t] = Butter, Cheby-1, Ellip, t=0,1,3');
disp('[-inf, 0, -20,  -1, 0, 5, 0, 4, t],  [-inf, 0, -20,  -1, 15, 20, 16, 20, t] = Cheby-2, t=2');
disp(' ');

P = input('enter design parameters (or, copy one of the above, 0 to exit) = ');

if P==0; return; end

switch length(P)          
   case 7,
      N = P(1); G0 = P(2); G = P(3); GB = P(4); f0 = P(5); Df = P(6); type = P(7); 
      w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;
      [w1,w2] = bandedge(w0,Dw);
      f1 = fs * w1/(2*pi); f2 = fs * w2/(2*pi);
   case 8,
      N = P(1); G0 = P(2); G = P(3); GB = P(4); f0 = P(5); Df = P(6); type = P(7); Gs = P(8); 
      w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;
      [w1,w2] = bandedge(w0,Dw);
      f1 = fs * w1/(2*pi); f2 = fs * w2/(2*pi);
   case 9,
      G0 = P(1); G = P(2); GB = P(3); Gs = P(4); 
      f1 = P(5); f2 = P(6); f1s = P(7); f2s = P(8); type = P(9);
      w1 = 2*pi*f1/fs; w2 = 2*pi*f2/fs; w1s = 2*pi*f1s/fs; w2s = 2*pi*f2s/fs;
      [w0,Dw] = bandedge(w1,w2,1); Dws = w2s - w1s;
      f0 = fs * w0/(2*pi);
      Nexact = hpeqord(G0,G,GB,Gs,Dw,Dws,type)
      N = ceil(Nexact);
end

r = mod(N,2);   

g0 = 10^(G0/20); g = 10^(G/20); gB = 10^(GB/20);        % gains in absolute units
if length(P)==8 | length(P)==9, gs = 10^(Gs/20); end

switch type,                                            % H0 and Hinf gains in absolute units
   case 0,
      h0 = g; hinf = g0;
   case 1,
      h0 = g^r * gB^(1-r); hinf = g0;
   case 2,
      h0 = g; hinf = g0^r * gB^(1-r);
   case 3,
      h0 = g^r * gB^(1-r); hinf = g0^r * gs^(1-r);
end

H0 = 20*log10(h0); Hinf = 20*log10(hinf);               % H0 and Hinf gains in dB     

if type==3,
   [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,3,Gs,tol)
else 
   [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,type)
end

f = linspace(0,20,1001); w = 2*pi*f/fs;

H = abs(fresp(B,A,w));

if (g0~=0) & (g~=0),             % equalizer + shelving filters

  H = 20*log10(H);

  if length(P)~=9,
     plot(f,H,'b-', [0,f1,f0,f2,20], [Hinf,GB,H0,GB,Hinf],'r.');  
  else
     plot(f,H,'b-', [0,f1s,f1,f0,f2,f2s,20], [Hinf,Gs,GB,H0,GB,Gs,Hinf],'r.');  
  end

  ylim([-15,15]); ytick(-15:3:15); ylabel('dB');
  xtick(0:2:20); xlabel('{\it f}  (kHz)'); grid;

else                              % ordinary bandpass, bandstop, lowpass, highpass 

  if length(P)~=9,
     plot(f,H,'b-', [0,f1,f0,f2,20], [hinf,gB,h0,gB,hinf], 'r.');  
  else
     plot(f,H,'b-', [0,f1s,f1,f0,f2,f2s,20], [hinf,gs,gB,h0,gB,gs,hinf],'r.');  
  end

  ylim([0,1.05]); ytick(0:0.1:1); ylabel('absolute units');
  xtick(0:2:20); xlabel('{\it f}  (kHz)'); grid;


end



