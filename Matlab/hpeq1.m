% hpeq1.m - high-order digital parametric equalizer design (using explicit design equations)
%
% Usage: [B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,type); 
%
%        [B,A,Bh,Ah] = hpeq1(N,G0,G,GB,w0,Dw);            Butterworth
%        [B,A,Bh,Ah] = hpeq1(N,G0,G,GB,w0,Dw,0);          Butterworth
%        [B,A,Bh,Ah] = hpeq1(N,G0,G,GB,w0,Dw,1);          Chebyshev-1
%        [B,A,Bh,Ah] = hpeq1(N,G0,G,GB,w0,Dw,2);          Chebyshev-2
%
% N  = analog filter order
% G0 = reference gain (all gains must be in dB, enter -Inf to get 0 in absolute units)
% G  = peak/cut gain
% GB = bandwidth gain
% w0 = peak/cut center frequency in units of radians/sample, i.e., w0=2*pi*f0/fs
% Dw = bandwidth in radians/sample, (if w0=pi, Dw = cutoff freq measured from Nyquist)
% type = 0,1,2, for Butterworth, Chebyshev-1, Chebyshev-2 (default is type=0)
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
%        computation of left and right bandedge frequencies:
%        for the passband: [w1,w2] = bandedge(w0,Dw)
%        for the stopband in the elliptic case: k = ellipdeg(N,k1,tol); WB = tan(Dw/2); 
%        Ws = WB/k; Dws = 2*atan(Ws), [w1s,w2s] = bandedge(w0,Dws)
%
%        produces same output as HPEQ, but it uses the explicit design equations

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

function [B,A,Bh,Ah] = hpeq1(N,G0,G,GB,w0,Dw,type)

if nargin==0, help hpeq1; return; end
if nargin==6, type=0; end

G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20); if type==3, Gs = 10^(Gs/20); end

r = rem(N,2); L = (N-r)/2;

Bh = [1 0 0]; Ah = [1 0 0]; B = [1 0 0 0 0]; A = [1 0 0 0 0];

if G==G0, return; end              % no filtering if G=G0

c0 = cos(w0); 

if w0==0,    c0=1;  end    % special cases
if w0==pi/2, c0=0;  end
if w0==pi,   c0=-1; end

WB = tan(Dw/2);
e = sqrt((G^2 - GB^2)/(GB^2 - G0^2)); 

g = G^(1/N); g0 = G0^(1/N); 

if type==0,                                         % Butterworth
  b = WB / e^(1/N);

  if r==1,                                       
    D = b + 1;
    Bh(1,:) = [g*b+g0, g*b-g0, 0]/D;                    
    Ah(1,:) = [1, (b-1)/D, 0];

    B(1,:) = [g0+g*b, -2*g0*c0, g0-g*b, 0, 0]/D;       
    A(1,:) = [1, [-2*c0, 1-b, 0, 0]/D];
  end

  for i=1:L,
    phi = (2*i-1)*pi/(2*N);                                
    si = sin(phi); 
    D = b^2 + 2*b*si + 1;
                           
    b0h = (g^2*b^2 + 2*g0*g*b*si + g0^2)/D;        
    b1h = 2*(g^2*b^2 - g0^2)/D;
    b2h = (g^2*b^2 - 2*g0*g*b*si + g0^2)/D;
    a1h = 2*(b^2 - 1)/D;
    a2h = (b^2 - 2*b*si + 1)/D;
    Bh(i+1,:) = [b0h, b1h, b2h];
    Ah(i+1,:) = [ 1,  a1h, a2h];
                                         
    b0 = (g^2*b^2 + g0^2 + 2*g*g0*si*b)/D;        
    b1 = -4*c0*(g0^2 + g*g0*si*b)/D;
    b2 = 2*((1+2*c0^2)*g0^2 - g^2*b^2)/D;
    b3 = -4*c0*(g0^2 - g*g0*si*b)/D;
    b4 = (g^2*b^2 + g0^2 - 2*g*g0*si*b)/D;
    a1 = -4*c0*(1 + si*b)/D;
    a2 = 2*(1+2*c0^2 - b^2)/D;
    a3 = -4*c0*(1 - si*b)/D;
    a4 = (b^2 - 2*si*b + 1)/D;
    B(i+1,:) = [b0, b1, b2, b3, b4];
    A(i+1,:) = [1,  a1, a2, a3, a4];
  end
end    % end type=0 

if type==1,                                         % Chebyshev-1
  ea = (1/e + sqrt(1+1/e^2))^(1/N);
  eb = (G/e + GB*sqrt(1+1/e^2))^(1/N);
  a = (ea - 1/ea)/2;                                % a=sinh(alpha), sinh(N*alpha)=1/e, ea=exp(alpha)			  
  b = (eb - g0^2/eb)/2;                             % b=g0*sinh(beta), sinh(N*beta)=G/(G0*e), eb=g0*exp(beta)

  if r==1,                                         
    D = a*WB + 1;
    Bh(1,:) = [b*WB+g0, b*WB-g0, 0]/D;          
    Ah(1,:) = [1, (a*WB-1)/D, 0];

    B(1,:) = [g0+b*WB, -2*g0*c0, g0-b*WB, 0, 0]/D;   
    A(1,:) = [1, [-2*c0, 1-a*WB, 0, 0]/D];
  end

  for i=1:L,
    phi = (2*i-1)*pi/(2*N);                                
    ci = cos(phi); si = sin(phi);
    D = WB^2*(a^2+ci^2) + 2*a*si*WB + 1;

    b0h = (WB^2*(b^2+g0^2*ci^2) + 2*g0*b*si*WB + g0^2)/D;     % 2nd order sections
    b1h = 2*(WB^2*(b^2+g0^2*ci^2) - g0^2)/D;
    b2h = (WB^2*(b^2+g0^2*ci^2) - 2*g0*b*si*WB + g0^2)/D;
    a1h = 2*(WB^2*(a^2+ci^2) - 1)/D;
    a2h = (WB^2*(a^2+ci^2) - 2*a*si*WB + 1)/D;
    Bh(i+1,:) = [b0h, b1h, b2h];
    Ah(i+1,:) = [ 1,  a1h, a2h];

    b0 = (WB^2*(b^2+g0^2*ci^2) + 2*g0*b*si*WB + g0^2)/D;      % 4th order sections 
    b1 = -4*c0*g0*(b*si*WB + g0)/D;
    b2 = 2*((1+2*c0^2)*g0^2 - WB^2*(b^2+g0^2*ci^2))/D;
    b3 = -4*c0*g0*(-b*si*WB + g0)/D;
    b4 = (WB^2*(b^2+g0^2*ci^2) - 2*g0*b*si*WB + g0^2)/D;
    a1 = -4*c0*(a*si*WB + 1)/D;
    a2 = 2*(1+2*c0^2 - WB^2*(a^2+ci^2))/D;
    a3 = -4*c0*(-a*si*WB + 1)/D;
    a4 = (WB^2*(a^2+ci^2) - 2*a*si*WB + 1)/D;
    B(i+1,:) = [b0, b1, b2, b3, b4];
    A(i+1,:) = [1,  a1, a2, a3, a4];
  end
end     % end type=1

if type==2,                                         % Chebyshev-2
  ea = (e + sqrt(e^2+1))^(1/N);
  eb = (G0*e + GB*sqrt(e^2+1))^(1/N); 
  a = (ea - 1/ea)/2;                                % a=sinh(alpha), sinh(N*alpha)=e; ea=exp(alpha)
  b = (eb - g^2/eb)/2;                              % b=g*sinh(beta), sinh(N*beta)=G0*e/G, eb=g*exp(beta)

  if r==1,                                       
    D = WB + a;
    Bh(1,:) = [g*WB+b, g*WB-b, 0]/D;            
    Ah(1,:) = [1, (WB-a)/D, 0];

    B(1,:) = [g*WB+b, -2*b*c0, b-g*WB, 0, 0]/D;     
    A(1,:) = [1, [-2*a*c0, a-WB, 0, 0]/D];
  end

  for i=1:L,
    phi = (2*i-1)*pi/(2*N);                                
    ci = cos(phi); si = sin(phi);
    D = WB^2 + 2*a*WB*si + a^2 + ci^2;
                                 
    b0h = (g^2*WB^2 + 2*WB*g*b*si + b^2 + g^2*ci^2)/D;   
    b1h = 2*(g^2*WB^2 - b^2 - g^2*ci^2)/D;
    b2h = (g^2*WB^2 - 2*WB*g*b*si + b^2 + g^2*ci^2)/D;
    a1h = 2*(WB^2 - a^2 - ci^2)/D;
    a2h = (WB^2 - 2*WB*a*si + a^2 + ci^2)/D;
    Bh(i+1,:) = [b0h, b1h, b2h];
    Ah(i+1,:) = [ 1,  a1h, a2h];
                                                
    b0 = (g^2*WB^2 + b^2 + g^2*ci^2 + 2*g*b*si*WB)/D;   
    b1 = -4*c0*(b^2 + g^2*ci^2 + g*b*si*WB)/D;
    b2 = 2*((1+2*c0^2)*(b^2+g^2*ci^2) - g^2*WB^2)/D;
    b3 = -4*c0*(b^2 + g^2*ci^2 - g*b*si*WB)/D;
    b4 = (g^2*WB^2 + b^2 + g^2*ci^2 - 2*g*b*si*WB)/D;
    a1 = -4*c0*(a^2+ci^2 + a*si*WB)/D;
    a2 = 2*((1+2*c0^2)*(a^2+ci^2) - WB^2)/D;
    a3 = -4*c0*(a^2 + ci^2 - a*si*WB)/D;
    a4 = (WB^2+a^2 + ci^2 - 2*a*si*WB)/D;
    B(i+1,:) = [b0, b1, b2, b3, b4];
    A(i+1,:) = [1,  a1, a2, a3, a4];
  end
end     % end type=2


if c0==1 | c0==-1 	        % LP or HP shelving filter
   B = Bh;                      % B,A are second-order
   A = Ah;
   B(:,2) = c0*B(:,2);	        % change sign if w0=pi
   A(:,2) = c0*A(:,2);
   B(:,4:5) = 0;                % make them (L+1)x5
   A(:,4:5) = 0;                % for convenience in using fresp
end  
