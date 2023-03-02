% decoupfilt.m - filtering in frequency-shifted 2nd-order cascaded normalized lattice sections
%
% Usage: [y,W] = decoupfilt(gamma,d,w0,x,W);
%            y = decoupfilt(gamma,d,w0,x);        (default start with W=0)
%
% gamma = Kx2 matrix whose rows are the reflection coefficients of lattice sections
% d     = Kx3 matrix whose rows are the ladder coefficients of decoupled sections
% w0    = center frequency (rads/sample)
% x     = vector of input samples
% W     = Kx8 matrix of input internal states
%
% y = vector of output samples
% W = Kx8 matrix of output internal states, after processing last x
%
% notes: may be used in conjunction with DIR2DECOUP to calculate the [gamma,d] coefficients
%
%        first section is always 1st order or 0th order, that is,
%        gamma(1,:) = [gamma(1,1), 0],  d(1,:) = [d(1,1), d(1,2), 0]
%
%        internal states are labeled as follows:
%           W(i,1:2) = contents of zhat delays for gamma1,tau1 section
%           W(i,3:4) = contents of zhat delays for gamma2,tau2 section
%           W(i,5:6) = contents of zhat delays for the leftmost (1+zhat^-1) factor
%           W(i,7:8) = contents of zhat delays for the rightmost (1+zhat^-1) factor

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

function [y,W] = decoupfilt(gamma,d,w0,x,W)

if nargin==0, help decoupfilt; return; end

K = size(gamma,1);         % number of sections

c0 = cos(w0); s0 = sqrt(1-c0^2);

if w0==0,    c0=1;  s0=0; end    % special cases
if w0==pi/2, c0=0;  s0=1; end
if w0==pi,   c0=-1; s0=0; end

tau = sqrt(1-gamma.^2);    % transmission coefficients

if nargin==4, W=zeros(K,8); end

y = x;

for n=1:length(x),
   v1 = W(1,2);                                  % first section, states W(1,:) = [W(1,1),W(1,2),0,0,0,0,0,0]
   x0 = tau(1,1)*y(n) - gamma(1,1)*v1;
   y1 = tau(1,1)*v1 + gamma(1,1)*y(n);
   y(n) = d(1,1)*y(n) + d(1,2)*y1;
   W(1,2) = c0*x0 - s0*W(1,1);
   W(1,1) = s0*x0 + c0*W(1,1);

   for i=2:K,                                       % there are 4 zhat delays, or, 8 internal states
      y1 = tau(i,2)*W(i,4) + gamma(i,2)*y(n);      
      x1 = tau(i,2)*y(n) - gamma(i,2)*W(i,4);
      v1 = tau(i,1)*W(i,2) + gamma(i,1)*x1;
      x0 = tau(i,1)*x1 - gamma(i,1)*W(i,2);
      v2 = x0 + W(i,6);                          
      y2 = v2 + W(i,8);
      y(n) = d(i,1)*y(n) + d(i,2)*y1 + d(i,3)*y2;
      W(i,4) = c0*v1 - s0*W(i,3);
      W(i,3) = s0*v1 + c0*W(i,3);
      W(i,2) = c0*x0 - s0*W(i,1);
      W(i,1) = s0*x0 + c0*W(i,1);
      W(i,8) = c0*v2 - s0*W(i,7);
      W(i,7) = s0*v2 + c0*W(i,7);
      W(i,6) = c0*x0 - s0*W(i,5);
      W(i,5) = s0*x0 + c0*W(i,5);
   end
end        




