% hpeq_a.m - high-order analog parametric equalizer design
%
% Usage: [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,type,Gs,tol); 
%
%        [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df);            Butterworth
%        [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,0);          Butterworth
%        [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,1);          Chebyshev-1
%        [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,2);          Chebyshev-2
%        [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,3,Gs,tol);   elliptic, e.g., tol = 1e-8
%        [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,3,Gs,M);     elliptic, tol = M = Landen iterations, e.g., M=5
%        [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,3,Gs);       elliptic, default tol = eps
%
% N  = analog filter order
% G0 = reference gain (all gains must be in dB, enter -Inf to get 0 in absolute units)
% G  = peak/cut gain 
% GB = bandwidth gain
% f0 = peak/cut center frequency in Hz, (f0=Inf, for highpass)
% Df = bandwidth in Hz,  (Df = cutoff frequency, for LP or HP)
% type = 0,1,2,3, for Butterworth, Chebyshev-1, Chebyshev-2, and elliptic (default is type=0)
% Gs = stopband gain, for elliptic case only
% tol = tolerance for elliptic case, e.g., tol = 1e-10, default value is tol = eps = 2.22e-16
%
% B,A   = rows are the numerator and denominator 4th order section coefficients of the equalizer
% Ba,Aa = rows are the numerator and denominator 2nd order coefficients of the lowpass shelving filter
%
% notes: G,GB,G0 are in dB and converted internally to absolute units, e.g. G => 10^(G/20).
%
%        gains must satisfy: G0<Gs<GB<G (boost), or, G0>Gs>GB>G (cut)  (exchange roles of Gs,GB for Cheby-2)
% 
%        B,A have size (L+1)x5, and Ba,Aa, size (L+1)x3, where L=floor(N/2)
%        when N is even, the first row is just a gain factor
%
%        computation of left and right bandedge frequencies:
%           [f1,f2] = bandedga(f0,Df);
%        for the stopband in the elliptic case:
%           k = ellipdeg(N,k1,tol); Dfs = Df/k; [f1s,f2s] = bandedga(f0,Dfs)
%
%        except for the definition of WB, Ba,Aa are identical to those of HPEQ
%
%        Run the scripts HPEQEX1_A, HPEQEX2_A to generate some examples
%
% see also, ELLIDPEG, ELLIPK, ASNE, CDE, BANDEDGE_A

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

function [B,A,Ba,Aa] = hpeq_a(N,G0,G,GB,f0,Df,type,Gs,tol)

if nargin==0, help hpeq_a; return; end
if nargin==6, type=0; end
if type==3 & nargin<=7, disp('must enter values for Gs,tol'); return; end
if type==3 & nargin==8, tol=eps; end

G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20); if type==3, Gs = 10^(Gs/20); end

r = rem(N,2); L = (N-r)/2;

w0 = 2*pi*f0; WB = 2*pi*Df;

if w0==Inf, WB = 1/WB; end    % LP->HP transformation is s = 1/p

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


switch w0     % frequency transformations: LP) s => s, HP) s => 1/s, and BP) s => s + w0^2/s
   case 0,                                
      B = Ba; 
      A = Aa;
   case Inf,                              
      B = fliplr(Ba);
      A = fliplr(Aa);
      if r==1,                                   % fix the first row
         B(1,:) = [fliplr(Ba(1,1:2)),0];
         A(1,:) = [fliplr(Aa(1,1:2)),0];
      else
         B(1,:) = [Ba(1,1), 0, 0];
         A(1,:) = [Aa(1,1), 0, 0];
      end
   otherwise,                              
      B = [w0^4*Ba(:,3), w0^2*Ba(:,2), Ba(:,1)+2*w0^2*Ba(:,3), Ba(:,2), Ba(:,3)];
      A = [w0^4*Aa(:,3), w0^2*Aa(:,2), Aa(:,1)+2*w0^2*Aa(:,3), Aa(:,2), Aa(:,3)];
      if r==1,
         B(1,:) = [w0^2*Ba(1,2), Ba(1,1), Ba(1,2), 0, 0];
         A(1,:) = [w0^2*Aa(1,2), Aa(1,1), Aa(1,2), 0, 0];
      else
         B(1,:) = [Ba(1,1), 0, 0, 0, 0];
         A(1,:) = [Aa(1,1), 0, 0, 0, 0];
      end
end






    











