% hpeqex1_a.m - example output from hpeq_a.m
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

tol = eps;

disp(' ');
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
disp('[5, 0, 12,    9, 6,4, 0],    [5, 0, 12,    9, 0,4, 0],    [5, 0, 12,    9, inf,14, 0]    = Butter  EQ,LP,HP');
disp('[5, 0, 12, 11.5, 6,4, 1],    [5, 0, 12, 11.5, 0,4, 1],    [5, 0, 12, 11.5, inf,14, 1]    = Cheby-1 EQ,LP,HP');
disp('[5, 0, 12,  0.5, 6,4, 2],    [5, 0, 12,  0.5, 0,4, 2],    [5, 0, 12,  0.5, inf,14, 2]    = Cheby-2 EQ,LP,HP');
disp('[5, 0, 12, 11.5, 6,4, 3, 1], [5, 0, 12, 11.5, 0,4, 3, 1], [5, 0, 12, 11.5, inf,14, 3, 1] = Ellip   EQ,LP,HP');
disp(' ');
disp('[0, 12, 11.5,  0.5, 6,8,5,9, 0], [0, 12, 11.5,  0.5, 0,4,0,5, 0], [0, 12, 11.5,  0.5, 14,inf,13,inf, 0] = Butter');
disp('[0, 12, 11.5,  0.5, 6,8,5,9, 1], [0, 12, 11.5,  0.5, 0,4,0,5, 1], [0, 12, 11.5,  0.5, 14,inf,13,inf, 1] = Cheby-1');
disp('[0, 12,  0.5, 11.5, 5,9,6,8, 2], [0, 12,  0.5, 11.5, 0,5,0,4, 2], [0, 12,  0.5, 11.5, 13,inf,14,inf, 2] = Cheby-2');
disp('[0, 12, 11.5,  0.5, 6,8,5,9, 3], [0, 12, 11.5,  0.5, 0,4,0,5, 3], [0, 12, 11.5,  0.5, 14,inf,13,inf, 3] = Ellip');

disp(' ');
disp('Ordinary BP, BS with fixed or unknown order:');
disp('-------------------------------------------- ');
disp('[5, -inf, 0,  -3, 6,4, 0],       [5, 0, -inf,  -3, 6,4, 0]     = Butter');
disp('[5, -inf, 0,  -1, 6,4, 1],       [5, 0, -inf, -20, 6,4, 1]     = Cheby-1');
disp('[5, -inf, 0, -20, 6,4, 2],       [5, 0, -inf,  -1, 6,4, 2]     = Cheby-2');
disp('[5, -inf, 0,  -1, 6,4, 3, -20],  [5, 0, -inf, -20, 6,4, 3, -1] = Ellip');
disp(' ');
disp('[-inf, 0,  -1, -20, 6,8,5,9, t],  [0, -inf, -20,  -1, 6,8,5,9, t] = Butter, Cheby-1, Ellip, t=0,1,3');
disp('[-inf, 0, -20,  -1, 5,9,6,8, t],  [0, -inf,  -1, -20, 5,9,6,8, t] = Cheby-2, t=2');

disp(' ');
disp('Ordinary LP, HP with fixed or unknown order:');
disp('-------------------------------------------- ');
disp('[5, -inf, 0,  -3, 0,4, 0],       [5, -inf, 0,  -3, inf,14, 0]      = Butter');
disp('[5, -inf, 0,  -1, 0,4, 1],       [5, -inf, 0,  -1, inf,14, 1]      = Cheby-1');
disp('[5, -inf, 0, -20, 0,4, 2],       [5, -inf, 0, -20, inf,14, 2]      = Cheby-2');
disp('[5, -inf, 0,  -1, 0,4, 3, -20],  [5, -inf, 0,  -1, inf,14, 3, -20] = Ellip');
disp(' ');
disp('[-inf, 0,  -1, -20, 0,4,0,5, t],  [-inf, 0,  -1, -20, 14,inf,13,inf, t] = Butter, Cheby-1, Ellip, t=0,1,3');
disp('[-inf, 0, -20,  -1, 0,5,0,4, t],  [-inf, 0, -20,  -1, 13,inf,14,inf, t] = Cheby-2, t=2');
disp(' ');

P = input('enter design parameters (or, copy one of the above, 0 to exit) = ');

if P==0; return; end

switch length(P)          
   case 7,
      N = P(1); G0 = P(2); G = P(3); GB = P(4); f0 = P(5); Df = P(6); type = P(7); 
      [f1,f2] = bandedge_a(f0,Df);
   case 8,
      N = P(1); G0 = P(2); G = P(3); GB = P(4); f0 = P(5); Df = P(6); type = P(7); Gs = P(8); 
      [f1,f2] = bandedge_a(f0,Df);
   case 9,
      G0 = P(1); G = P(2); GB = P(3); Gs = P(4); 
      f1 = P(5); f2 = P(6); f1s = P(7); f2s = P(8); type = P(9);
      [f0,Df] = bandedge_a(f1,f2,1); [f0s, Dfs] = bandedge_a(f1s,f2s,1);
      if f0==inf,                
         Nexact = hpeqord_a(G0,G,GB,Gs,1/Df,1/Dfs,type)    % highpass case, Df=f1, Dfs=f1s, f2=f2s=Inf
      else
         Nexact = hpeqord_a(G0,G,GB,Gs,Df,Dfs,type)
      end
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
   [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,3,Gs,tol)
else 
   [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,type)
end

f = linspace(0,20,1001); 

H = abs(fresp_a(B,A,f));

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



