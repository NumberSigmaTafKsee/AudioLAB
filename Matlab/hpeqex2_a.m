% hpeqex2_a.m - example output from hpeq_a.m
%
% high-order designs with given 3-dB width
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

tol = eps;

disp(' ');
disp('Equalizer designs based on 3-dB width')
disp('-------------------------------------');
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
disp('[5, 0, 12, 9,    9, 6,4, 0],      [5, 0, 12, 9,    9, 0,4, 0],      [5, 0, 12, 9,    9, inf,14, 0]    = Butter  EQ,LP,HP');
disp('[5, 0, 12, 11.5, 9, 6,4, 1],      [5, 0, 12, 11.5, 9, 0,4, 1],      [5, 0, 12, 11.5, 9, inf,14, 1]    = Cheby-1 EQ,LP,HP');
disp('[5, 0, 12, 0.5,  9, 6,4, 2],      [5, 0, 12, 0.5,  9, 0,4, 2],      [5, 0, 12, 0.5,  9, inf,14, 2]    = Cheby-2 EQ,LP,HP');
disp('[5, 0, 12, 11.5, 9, 6,4, 3, 0.5], [5, 0, 12, 11.5, 9, 0,4, 3, 0.5], [5, 0, 12, 11.5, 9, inf,14, 3, 1] = Ellip   EQ,LP,HP');

disp(' ');
disp('EQ examples with unknown order (set f1=f1b=0 for LP, or, f2=f2b=Inf for HP):');
disp('----------------------------------------------------------------------------');
disp('[1, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 0],      [2, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 0]      = Butter,   m=1,2');
disp('[1, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 1],      [2, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 1]      = Cheby-1,  m=1,2');
disp('[1, 0, 12, 0.1,  9, 5.5, 8.5, 6,   8,   2],      [2, 0, 12, 0.1,  9, 5.5, 8.5, 6,   8,   2]      = Cheby-2,  m=1,2');
disp('[1, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 3, 0.1], [2, 0, 12, 11.9, 9, 6,   8,   5.5, 8.5, 3, 0.1] = Ellip,    m=1,2');
disp('[1, 0, 12, 11.9, 9, 13, inf, 12, inf, 3, 0.1],   [2, 0, 12, 11.9, 9, 13, inf, 12, inf, 3, 0.1]   = Ellip HP, m=1,2');

disp(' ');
disp('Ordinary BP, BS with fixed order (set f0=0,Inf for LP,HP):');
disp('---------------------------------------------------------- ');
disp('[5, -inf, 0, -1,  -3, 6, 4, 0],      [5, 0, -inf, -20, -3, 6, 4, 0]     = Butter  EQ,LP,HP');
disp('[5, -inf, 0, -1,  -3, 6, 4, 1],      [5, 0, -inf, -20, -3, 6, 4, 1]     = Cheby-1 EQ,LP,HP');
disp('[5, -inf, 0, -20, -3, 6, 4, 2],      [5, 0, -inf, -1,  -3, 6, 4, 2]     = Cheby-2 EQ,LP,HP');
disp('[5, -inf, 0, -1,  -3, 6, 4, 3, -20], [5, 0, -inf, -20, -3, 6, 4, 3, -1] = Ellip   EQ,LP,HP');

disp(' ');
disp('Ordinary BP, BS with unknown order (set f1=f1b=0 for LP, or, f2=f2b=Inf for HP):');
disp('-------------------------------------------------------------------------------- ');
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
   case 9,
      N = P(1); G0 = P(2); G = P(3); GB = P(4); Gb = P(5); f0 = P(6); Dfb = P(7); type = P(8); Gs = P(9);
   case 10,
      m = P(1); G0 = P(2); G = P(3); GB = P(4); Gb = P(5); 
      f1 = P(6); f2 = P(7); f1b = P(8); f2b = P(9); type = P(10);
      [f0,Df] = bandedge_a(f1,f2,1);  [f0b,Dfb] = bandedge_a(f1b,f2b,1);
   case 11,
      m = P(1); G0 = P(2); G = P(3); GB = P(4); Gb = P(5); 
      f1 = P(6); f2 = P(7); f1b = P(8); f2b = P(9); type = P(10); Gs = P(11);
      [f0,Df] = bandedge_a(f1,f2,1);  [f0b,Dfb] = bandedge_a(f1b,f2b,1);
end

switch length(P)
  case 8,                                                % type=0,1,2
    if f0==Inf, Dfb = 1/Dfb; end                         % 1/s transformation for HP case 
    Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,type);                 % Df = bandwidth at GB
    if f0==inf, Df=1/Df; Dfb = 1/Dfb; end                % correct HP case
    [f1,f2] = bandedge_a(f0,Df);                           % bandedges at GB
    [f1b,f2b] = bandedge_a(f0,Dfb);                        % bandedges at Gb
  case 9,                                                % type=3
    if f0==Inf, Dfb = 1/Dfb; end      
    Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,3,Gs);
    if f0==inf, Df=1/Df; Dfb = 1/Dfb; end
    [f1,f2] = bandedge_a(f0,Df); 
    [f1b,f2b] = bandedge_a(f0,Dfb);
  case 10,                                               % type=0,1,2
    if f0==inf, Df=1/Df; Dfb = 1/Dfb; end
    N = ceil(hpeqord_a(G0,G,GB,Gb,Df,Dfb,type));          % filter order from Df,Dfb
    if m==2, 
       Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,type); f0 = f0b;    % match Dfb, otherwise, match Df
    end
    if f0==inf, Df=1/Df; Dfb = 1/Dfb; end
  case 11,                                               % type=3
    if f0==inf, Df=1/Df; Dfb = 1/Dfb; end
    N = ceil(hpeqord_a(G0,G,GB,Gs,Df,Dfb,3,Gb));
    if m==2, 
       Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,3,Gs); f0 = f0b; 
    end
    if f0==inf, Df=1/Df; Dfb = 1/Dfb; end
end

g0 = 10^(G0/20); g = 10^(G/20); gB = 10^(GB/20); gb = 10^(Gb/20);    
if length(P)==9 | length(P)==11, gs = 10^(Gs/20); end

r = mod(N,2);   

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
   [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,3,Gs,tol)    
else 
   [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,type)
end

figure;

f = linspace(0,20,1001); 

H = abs(fresp_a(B,A,f));

if (g0~=0) & (g~=0),         % equalizer + shelving filters, y-scale in dB

  H = 20*log10(H);
  plot(f,H,'b-', [0,f1b,f1,f0,f2,f2b,20], [Hinf,Gb,GB,H0,GB,Gb,Hinf],'r.');  
  ylim([-15,15]); ytick(-15:3:15); ylabel('dB');
  xtick(0:2:20); xlabel('{\it f}  (kHz)'); grid;

else                         % ordinary bandpass, bandstop, lowpass, highpass 

  plot(f,H,'b-', [0,f1b,f1,f0,f2,f2b,20], [hinf,gb,gB,h0,gB,gb,hinf], 'r.');  
  ylim([0,1.05]); ytick(0:0.1:1); ylabel('absolute units');
  xtick(0:2:20); xlabel('{\it f}  (kHz)'); grid;

end

