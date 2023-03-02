% hpeq.m - high-order digital parametric equalizer design
%
% Usage: [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,type,Gs,tol); 
%
%        [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw);            Butterworth
%        [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,0);          Butterworth
%        [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,1);          Chebyshev-1
%        [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,2);          Chebyshev-2
%        [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,3,Gs,tol);   elliptic, e.g., tol = 1e-8
%        [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,3,Gs,M);     elliptic, tol = M = Landen iterations, e.g., M=5
%        [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,3,Gs);       elliptic, default tol = eps
%
% N  = analog filter order
% G0 = reference gain (all gains must be in dB, enter -Inf to get 0 in absolute units)
% G  = peak/cut gain
% GB = bandwidth gain
% w0 = peak/cut center frequency in units of radians/sample, i.e., w0=2*pi*f0/fs
% Dw = bandwidth in radians/sample, (if w0=pi, Dw = cutoff freq measured from Nyquist)
% type = 0,1,2,3, for Butterworth, Chebyshev-1, Chebyshev-2, and elliptic (default is type=0)
% Gs = stopband gain, for elliptic case only
% tol = tolerance for elliptic case, e.g., tol = 1e-10, default value is tol = eps = 2.22e-16
%
% B,A   = rows are the numerator and denominator 4th order section coefficients of the equalizer
% Bh,Ah = rows are the numerator and denominator 2nd order coefficients of the lowpass shelving filter
%
% notes: G,GB,G0 are in dB and converted internally to absolute units, e.g. G => 10^(G/20)
%
%        gains must satisfy: G0<Gs<GB<G (boost), or, G0>Gs>GB>G (cut)  (exchange roles of Gs,GB for Cheby-2)
% 
%        w0 = 2*pi*f0/fs, Dw = 2*pi*Df/fs, with f0,Df,fs in Hz
% 
%        B,A have size (L+1)x5, and Bh,Ah, size (L+1)x3, where L=floor(N/2)
%        when N is even, the first row is just a gain factor
%
%        left and right bandedge frequencies: [w1,w2] = bandedge(w0,Dw)
%        for the stopband in the elliptic case: [w1s,w2s] = bandedge(w0,Dws), where 
%        k = ellipdeg(N,k1,tol); WB = tan(Dw/2); Ws = WB/k; Dws = 2*atan(Ws)
%
%        Run the functions HPEQEX0, HPEQEX1, HPEQEQ2 to generate some examples
%
% see also, BLT, ELLIDPEG, ELLIPK, ASNE, CDE, BANDEDGE, HPEQBW, OCTBW

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


function [B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type, Gs, tol)

if nargin==0, help hpeq; return; end
if nargin==6, type=0; end
if type==3 & nargin<=7, disp('must enter values for Gs,tol'); return; end
if type==3 & nargin==8, tol=eps; end

G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20); if type==3, Gs = 10^(Gs/20); end

r = rem(N,2); L = (N-r)/2;

%Bh = [1 0 0]; Ah = [1 0 0]; A = [1 0 0 0 0]; B = [1 0 0 0 0];

if G==G0,                             % no filtering if G=G0
   Bh = G0*[1 0 0];     Ah = [1 0 0]; 
   B  = G0*[1 0 0 0 0]; A  = [1 0 0 0 0];
   return; 
end              

c0 = cos(w0); 

if w0==0,    c0=1;  end    % special cases
if w0==pi/2, c0=0;  end
if w0==pi,   c0=-1; end

WB = tan(Dw/2);
e = sqrt((G^2 - GB^2)/(GB^2 - G0^2)); 

g = G^(1/N); g0 = G0^(1/N); 

switch type
  case 0,
    a = e^(1/N);
    b = g0*a;
  case 1,
    eu = (1/e + sqrt(1+1/e^2))^(1/N);
    ew = (G/e + GB*sqrt(1+1/e^2))^(1/N);
    a = (eu - 1/eu)/2;			
    b = (ew - g0^2/ew)/2;	              
  case 2,
    eu = (e + sqrt(1+e^2))^(1/N);
    ew = (G0*e + GB*sqrt(1+e^2))^(1/N);
    a = (eu - 1/eu)/2;
    b = (ew - g^2/ew)/2;
  case 3,
    es = sqrt((G^2 - Gs^2)/(Gs^2 - G0^2)); 
    k1 = e/es;
    k = ellipdeg(N, k1, tol);
    if G0~=0, 
       ju0 = asne(j*G/e/G0, k1, tol)/N;    % not used when G0=0
    end    
    jv0 = asne(j/e, k1, tol)/N;
end

if r==0,
  switch type
    case {0,1,2}
      Ba(1,:) = [1, 0, 0]; 
      Aa(1,:) = [1, 0, 0];
    case 3
      Ba(1,:) = [1, 0, 0] * GB;
      Aa(1,:) = [1, 0, 0];
   end
end 
     
if r==1,
  switch type
    case 0
      Ba(1,:) = [g*WB, b, 0]; 
      Aa(1,:) = [WB,   a, 0];
    case 1
      Ba(1,:) = [b*WB, g0, 0]; 
      Aa(1,:) = [a*WB, 1,  0]; 
    case 2
      Ba(1,:) = [g*WB, b, 0]; 
      Aa(1,:) = [WB,   a, 0];
    case 3
      if G0==0 & G~=0,
         B00 = G*WB; B01 = 0;
      elseif G0~=0 & G==0,
         K=ellipk(k,tol); K1=ellipk(k1,tol); 
         B00 = 0; B01 = G0*e*N*K1/K;
      else                                  % G0~=0 and G~=0
         z0 = real(j*cde(-1+ju0,k,tol));    % it's supposed to be real
         B00 = G*WB; B01 = -G/z0;
      end    
      p0 = real(j*cde(-1+jv0,k,tol));
      A00 = WB; A01 = -1/p0;
      Ba(1,:) = [B00,B01,0];
      Aa(1,:) = [A00,A01,0];
   end
end 

if L>0, 
   i = (1:L)';
   ui = (2*i-1)/N;  
   ci = cos(pi*ui/2); si = sin(pi*ui/2);
   v = ones(L,1);

   switch type
      case 0,
        Ba(1+i,:) = [g^2*WB^2*v, 2*g*b*si*WB, b^2*v];
        Aa(1+i,:) = [WB^2*v, 2*a*si*WB, a^2*v];
      case 1,
        Ba(1+i,:) = [WB^2*(b^2+g0^2*ci.^2), 2*g0*b*si*WB, g0^2*v];
        Aa(1+i,:) = [WB^2*(a^2+ci.^2), 2*a*si*WB, v];
      case 2,
        Ba(1+i,:) = [g^2*WB^2*v, 2*g*b*si*WB, b^2+g^2*ci.^2];
        Aa(1+i,:) = [WB^2*v, 2*a*si*WB, a^2+ci.^2];
      case 3,
        if G0==0 & G~=0,
           zeros = j ./ (k*cde(ui,k,tol));
        elseif G0~=0 & G==0,
           zeros = j*cde(ui,k,tol);
        else                               % G0~=0 and G~=0
           zeros = j*cde(ui-ju0,k,tol);
        end
        poles = j*cde(ui-jv0,k,tol);
        Bi0 = WB^2*v; Bi1 = -2*WB*real(1./zeros); Bi2 = abs(1./zeros).^2; 
        Ai0 = WB^2*v; Ai1 = -2*WB*real(1./poles); Ai2 = abs(1./poles).^2;
        Ba(1+i,:) = [Bi0, Bi1, Bi2];
        Aa(1+i,:) = [Ai0, Ai1, Ai2];
   end
end

[B,A,Bh,Ah] = blt(Ba,Aa,w0);







    











