% bandwidth.m - calculate bandwidth at any level
%
% Usage: Dwb = bandwidth(N,G0,G,GB,Gb,Dw,type,Gs); 
%
%        Dwb = bandwidth(N,G0,G,GB,Gb,Dw);          Butterworth (equivalent to type=0)
%        Dwb = bandwidth(N,G0,G,GB,Gb,Dw,0);        Butterworth
%        Dwb = bandwidth(N,G0,G,GB,Gb,Dw,1);        Chebyshev-1
%        Dwb = bandwidth(N,G0,G,GB,Gb,Dw,2);        Chebyshev-2
%        Dwb = bandwidth(N,G0,G,GB,Gb,Dw,3,Gs);     Elliptic
%
% N  = analog filter order
% G0 = reference gain (all gains must be in dB, enter -inf to get 0 in absolute units)
% G  = peak/cut gain 
% GB = bandwidth gain
% Gb = intermediate gain level (e.g., 3-dB below G) at which bandwidth is to be calculated
% Dw = design bandwidth at level GB (in units of rads/sample)
% type = 0,1,2,3, for Butterworth, Chebyshev-1, Chebyshev-2, and Elliptic (default is type=0)
% Gs = stopband gain in dB, for elliptic case only
%
% Dwb = bandwidth at level Gb (in rads/sample)
%
% notes: it solves the magnitude equation: (G^2 + G0^2*e^2*FN(wb)^2)/(1 + e^2*FN(wb)^2) = Gb^2,
%        which is equivalent to FN(wb) = eb/e, for wb = Wb/WB, where Wb = tan(Dwb/2) and WB = tan(Dw/2), 
%        and then computes Dwb = 2*atan(WB*wb)
%
%        it is essentially the reverse of HPEQBW, which calculates Dw from Dwb
%
%        special cases: Gb = G - 3*sign(G)  ==>  Dwb = 3-dB width (3-dB below peak or above cut)
%                       Gb = Gs  ==>  Dwb = Dws = width at stopband level
%
%        boost: G0<Gs<Gb<GB<G  (type=0,1,3),  G0<GB<Gb<G (type=2)
%        cut:   G0>Gs>Gb>GB>G  (type=0,1,3),  G0>GB>Gb>G (type=2)
%
%        to calculate the bandedge frequencies from Dwb, use BANDEDGE: [w1b,w2b] = bandedge(w0,Dwb)

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

function Dwb = bandwidth(N,G0,G,GB,Gb,Dw,type,Gs)

if nargin==0, help bandwidth; return; end
if nargin==6, type=0; end
if type==3 & nargin==7, disp('must enter a value for Gs'); return; end

G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20); Gb = 10^(Gb/20);
if nargin==8, Gs = 10^(Gs/20); end

e = sqrt((G^2 - GB^2)/(GB^2 - G0^2)); 
eb = sqrt((G^2 - Gb^2)/(Gb^2 - G0^2));
 
Fb = eb/e;                % Fb = FN(wb), where wb = Wb/WB

switch type,              % solve equation FN(wb) = eb/e, for wb
   case 0,
      wb = Fb^(1/N);
   case 1,
      u = acos(Fb)/N; wb = cos(u);
   case 2,
      u = acos(1/Fb)/N; wb = 1/cos(u);
   case 3,
      tol = eps;                              % may be changed, e.g., tol=1e-15, or, tol=5 Landen iterations
      es = sqrt((G^2 - Gs^2)/(Gs^2 - G0^2)); 
      k1 = e/es; 
      k = ellipdeg(N,k1,tol);
      u = acde(Fb,k1,tol)/N; wb = cde(u,k,tol);   % solves cd(N*u*K1,k1) = Fb for u, and then, wb = cd(u*K,k)
end

WB = tan(Dw/2);
Wb = WB*wb;  
Dwb = 2*atan(Wb);












