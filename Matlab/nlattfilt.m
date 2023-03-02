% nlattfilt.m - filtering in frequency-shifted 2nd-order cascaded normalized lattice sections
%
% Usage: [y,W] = nlattfilt(gamma,d,w0,x,W);
%            y = nlattfilt(gamma,d,w0,x);        (default start with W=0)
%
% gamma = Kx2 matrix whose rows are the reflection coefficients of lattice sections
% d     = Kx3 matrix whose rows are the ladder coefficients of lattice sections
% w0    = center frequency (rads/sample)
% x     = vector of input samples
% W     = Kx4 matrix of input internal states
%
% y = vector of output samples
% W = Kx4 matrix of output internal states, after processing last x
%
% notes: there are K lattice sections with 4 delays each
%        for the i-th section, the delay contents are [W(i,1),W(i,2),W(i,3),W(i,4]]
%        labeled such that W(i,4) is closest to the output
%
%        may be used in conjunction with DIR2LATT to calculate the [gamma,d] coefficients

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

function [y,W] = nlattfilt(gamma,d,w0,x,W)

if nargin==0, help nlattfilt; return; end

K = size(gamma,1);         % number of sections

c0 = cos(w0); s0 = sqrt(1-c0^2);

if w0==0,    c0=1;  s0=0; end    % special cases
if w0==pi/2, c0=0;  s0=1; end
if w0==pi,   c0=-1; s0=0; end

tau = sqrt(1-gamma.^2);    % transmission coefficients

if nargin==4, W=zeros(K,4); end

y = x;

for n=1:length(x),                                     
   for i=1:K,                                          
      y2 = gamma(i,2) * y(n) + tau(i,2) * W(i,4);        
      x1 = tau(i,2) * y(n) - gamma(i,2) * W(i,4);
      y1 = gamma(i,1) * x1 + tau(i,1) * W(i,2);
      y0 = tau(i,1) * x1 - gamma(i,1) * W(i,2);
      y(n) = d(i,1) * y0 + d(i,2) * y1 + d(i,3) * y2;
      W(i,2) = c0 * y0 - s0 * W(i,1);
      W(i,1) = s0 * y0 + c0 * W(i,1);
      W(i,4) = c0 * y1 - s0 * W(i,3);
      W(i,3) = s0 * y1 + c0 * W(i,3);
   end
end
   


