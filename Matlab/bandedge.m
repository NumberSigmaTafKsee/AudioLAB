% bandedge.m - calculate left and right bandedge frequencies from bilinear transformation
%
% Usage: [w1,w2] = bandedge(w0,Dw,dir);
%    
%        [w1,w2] = bandedge(w0,Dw);          equivalent to dir=0
%        [w1,w2] = bandedge(w0,Dw,0);        calculate w1,w2 from w0,Dw
%        [w0,Dw] = bandedge(w1,w2,1);        calculate w0,Dw from w1,w2
%
% w0  = center frequency in radians/sample (w0 = 2*pi*f0/fs)
% Dw  = bandwidth in radians per sample
% dir = 0,1, direction of calculation
%
% w1,w2 = left and right bandedge frequencies in rads/sample 
%
% notes: dir=0 case computes w1,w2 from w0,Dw as follows:
%            WB = tan(Dw/2); c0=cos(w0); s0=sin(w0);
%            cos(w1) = (c0 + WB*sqrt(WB^2+s0^2))/(WB^2+1);
%            cos(w2) = (c0 - WB*sqrt(WB^2+s0^2))/(WB^2+1);
%
%        dir=1 case computes w0,Dw from w1,w2 as follows:
%            Dw = w2-w1; cos(w0) = sin(w1+w2)/(sin(w1)+sin(w2));
%
%        LP case: w0=0,  Dw=cutoff measured from DC,      w1=0,     w2=Dw
%        HP case: w0=pi, Dw=cutoff measured from Nyquist, w1=pi-Dw, w2=pi

% --------------------------------------------------------------
% Copyright (c) 2005 by Sophocles J. Orfanidis
% 
% Address: Sophocles J. Orfanidis                       
%          ECE Department, Rutgers University          
%          94 Brett Road, Piscataway, NJ 08854-8058
% Email:   orfanidi@ece.rutgers.edu
% Date:    June 15, 2005 (revised June 19, 2006)
% 
% Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
%
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% 
% tested with MATLAB R11.1 and R14
% --------------------------------------------------------------

function [w1,w2] = bandedge(w0,Dw,dir);

if nargin==0, help bandedge; return; end
if nargin==2, dir=0; end

if dir==0,
   WB = tan(Dw/2); c0 = cos(w0); s0=sin(w0);
   w1 = acos((c0 + WB*sqrt(WB^2+s0^2))/(1+WB^2));
   w2 = acos((c0 - WB*sqrt(WB^2+s0^2))/(1+WB^2));
else
   w2 = Dw-w0;
   if Dw==pi, 
      w1=pi;
   else 
      w1 = acos(sin(w0+Dw)/(sin(w0)+sin(Dw)));
   end
end






