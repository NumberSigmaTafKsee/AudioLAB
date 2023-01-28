% hpeqord.m - order determination of digital parametric equalizer
%
% Usage: [N,k,k1] = hpeqord(G0,G,GB,Gs,Dw,Dws,type,Gb); 
%
%        [N,k,k1] = hpeqord(G0,G,GB,Gs,Dw,Dws);          Butterworth
%        [N,k,k1] = hpeqord(G0,G,GB,Gs,Dw,Dws,0);        Butterworth
%        [N,k,k1] = hpeqord(G0,G,GB,Gs,Dw,Dws,1);        Chebyshev-1
%        [N,k,k1] = hpeqord(G0,G,GB,Gs,Dw,Dws,2);        Chebyshev-2
%        [N,k,k1] = hpeqord(G0,G,GB,Gs,Dw,Dws,3);        elliptic
%        [N,k,k1] = hpeqord(G0,G,GB,Gs,Dw,Dwb,3,Gb);     elliptic, Dwb at Gb
%
% G0  = reference gain (all gains must be in dB, enter -Inf to get 0 in absolute units)
% G   = peak/cut gain
% GB  = passband bandwidth gain
% Gs  = stopband gain
% Dw  = bandwidth in [radians/sample] at level GB
% Dws = bandwidth in [radians/sample] at level Gs (or, at Gb)
% type = 0,1,2,3, for Butterworth, Chebyshev-1, Chebyshev-2, and elliptic (default is type=0)
% Gb  = intermediate bandwidth level, elliptic case only, with Dwb measured at Gb not Gs
%
% N    = exact solution for the analog filter order, not rounded-up to the next integer
% k,k1 = filter order design parameters 
%
% notes: solves the degree equation for the filter order N
%        N is the exact non-integer solution and must be rounded up (or down, is so desired)
%
%        for Butterworth, Chebyshev-1, and elliptic:
%              G>GB>Gs>G0 and Dw<Dws, for a boost   (elliptic Gb case: G>GB>Gb>Gs>G0)
%              G<GB<Gs<G0 and Dw<Dws, for a cut     
%        for Chebyshev-2:
%              G>Gs>GB>G0 and Dw>Dws, for a boost
%              G<Gs<GB<G0 and Dw>Dws, for a cut
%
%        given the bandedge frequencies [w1,w2] and [w1s,w2s], use the following steps to design the filter:
%
%           Dw = w2-w1; Dws = w2s-w1s;
%           N = hpeqord(G0,G,GB,Gs,Dw,Dws,type,Gb); N = ceil(N);
%           w0 = bandedge(w1,w2,1); 
%           [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,type,Gs,tol);

% --------------------------------------------------------------
% Copyright (c) 2005 by Sophocles J. Orfanidis
% 
% Address: Sophocles J. Orfanidis                       
%          ECE Department, Rutgers University          
%          94 Brett Road, Piscataway, NJ 08854-8058
% Email:   orfanidi@ece.rutgers.edu
% Date:    June 15, 2005
% 
% Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
%
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% 
% tested with MATLAB R11.1 and R14
% --------------------------------------------------------------

function [N,k,k1] = hpeqord(G0, G, GB, Gs, Dw, Dws, type, Gb)

if nargin==0, help hpeqord; return; end
if nargin==6, type=0; end
if nargin==8, Nmax = 15; end    % max filter order to test

G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20); Gs = 10^(Gs/20);
if nargin==8, Gb = 10^(Gb/20); end

WB = tan(Dw/2); 
Ws = tan(Dws/2);
k = WB/Ws;

e = sqrt((G^2 - GB^2)/(GB^2 - G0^2)); 
es = sqrt((G^2 - Gs^2)/(Gs^2 - G0^2)); 
k1 = e/es;

switch type
  case 0,
     N = log(k1)/log(k);
  case 1,
     N = acosh(1/k1)/acosh(1/k);
  case 2,
     N = acosh(k1)/acosh(k);
  case 3,
     tol = eps;                  % may be changed, e.g., tol=1e-15, or, tol=5 Landen iterations
     if nargin==7,                                   % N based on GB,Gs
        [K,Kprime] = ellipk(k,tol);
        [K1,K1prime] = ellipk(k1,tol);
        N = (K1prime/K1)/(Kprime/K);
     else                                            % Gb-case
        eb = sqrt((G^2 - Gb^2)/(Gb^2 - G0^2));
        Fb = eb/e;                                   % Fb is F_N(wb)
        wb = 1/k;                                    % wb = tan(Dws/2)/tan(Dw/2), Dws plays role of Dwb
        for N=1:Nmax,                                % keep testing successive N's
           k = ellipdeg(N,k1,tol);                
           u = acde(Fb,k1,tol)/N;                    % solve cd(N*u*K1,k1) = eb/e for u
           err = cde(u,k) - wb;                      % error keeps decreasing from positive values
           if err<=0, break; end                     % N is lowest such that err<=0
        end      
     end
     if N==1,                                                   % N=1 because wb>eb/e 
        Gb2 = (G^2 + G0^2*e^2*wb^2)/(1 + e^2*wb^2);
        fprintf('\n\nhpeqord: \n');
        fprintf('   condition wb < eb/e failed, because wb = %6.4f,  eb/e = %6.4f\n', wb,Fb);
        fprintf('   Gb = %4.2f dB was too close to G = %4.2f dB\n', 20*log10([Gb,G]));
        fprintf('   G,Gb must be separated by |G-Gb| > %6.4f dB \n', 10*log10(G^2/Gb2));
        fprintf('   N was set to N = 1, but higher values of N will work, e.g., try N = 2-5\n\n');
     end
end


