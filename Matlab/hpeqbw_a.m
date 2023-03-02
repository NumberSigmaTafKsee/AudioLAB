% hpeqbw_a.m - remap bandwidth of high-order analog parametric equalizer
%
% Usage: Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,type,Gs,tol); 
%
%        Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb);          Butterworth (equivalent to type=0)
%        Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,0);        Butterworth
%        Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,1);        Chebyshev-1
%        Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,2);        Chebyshev-2
%        Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,3,Gs);     Elliptic
%
% N   = analog filter order
% G0  = reference gain (all gains must be in dB, enter -Inf to get 0 in absolute units)
% G   = peak/cut gain
% GB  = bandwidth gain
% Gb  = intermediate gain (e.g., 3-dB below G)
% Dfb = bandwidth at level Gb in Hz
% type = 0,1,2,3, for Butterworth, Chebyshev-1, Chebyshev-2, and elliptic (default is type=0)
% Gs = stopband gain, for elliptic case only
%
% Df = bandwidth at level GB in Hz
%
% notes: given bandwidth Dfb at level Gb, it computes Df at the design level GB
%
%        Df may be used in HPEQA to design the filter, that is, 
%        [B,A,Bh,Ah] = hpeqa(N,G0,G,GB,f0,Df,type,Gs,tol);
%
%        it solves the magnitude equation: (G^2 + G0^2*e^2*FN(wb)^2)/(1 + e^2*FN(wb)^2) = Gb^2  
%        for wb = Dfb/Df, and computes Df = 2*atan(Dfb/wb)
%        
%        boost: G0<Gs<Gb<GB<G  (type=0,1,3),  G0<GB<Gb<G (type=2)
%        cut:   G0>Gs>Gb>GB>G  (type=0,1,3),  G0>GB>Gb>G (type=2)
%
%        example: G0=0; G=12; GB=11.99; Gs=0.01; Gb=9 = 3-dB below peak   
%                 N=4; type = 3; f0 = 4; Dfb = 3; [f1,f2] = bandedge_a(f0,Dfb);
%                 Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,3,Gs); [B,A,Bh,Ah] = hpeq_a(N,G0,G,GB,f0,Df,3,Gs);
%                 f=linspace(0,10,1001); H=20*log10(abs(fresp_a(B,A,f))); 
%                 plot(f,H, [f1,f2],[Gb,Gb],'r.'); grid; ytick(0:1:12); xtick(0:1:10);
%
%        it uses the functions ELLIPDEG, ACDE, CDE

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
% References: R. A. Losada and V. Pellissier, "Designing IIR Filters with a Given 
%             3-dB Point," IEEE Signal Processing Mag., vol.22, no.4, 95, July 2005.
%
%             Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%             Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
%
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% 
% tested with MATLAB R11.1 and R14
% -------------------------------------------------------------------------

function Df = hpeqbw_a(N,G0,G,GB,Gb,Dfb,type,Gs,tol)

if nargin==0, help hpeqbw_a; return; end
if nargin==6, type=0; end
if type==3 & nargin==7, disp('must enter value for Gs'); return; end

G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20); Gb = 10^(Gb/20);
if nargin==8, Gs = 10^(Gs/20); end

e = sqrt((G^2-GB^2)/(GB^2-G0^2)); 
eb = sqrt((G^2-Gb^2)/(Gb^2-G0^2));
 
Fb = eb/e;                % Fb = FN(wb), where wb = Dfb/Df

switch type,
   case 0,
      wb = Fb^(1/N);
   case 1,
      u = acos(Fb)/N; 
      wb = cos(u);
   case 2,
      u = acos(1/Fb)/N; 
      wb = 1/cos(u);
   case 3,
      tol = eps;                              % may be changed, e.g., tol=1e-15, or, tol=5 Landen iterations
      es = sqrt((G^2-Gs^2)/(Gs^2-G0^2)); 
      k1 = e/es; 
      k = ellipdeg(N,k1,tol);
      u = acde(Fb,k1,tol)/N; 
      wb = cde(u,k,tol);
end

Df = Dfb/wb;  













