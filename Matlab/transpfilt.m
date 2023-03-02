% transpfilt.m - filtering in frequency-transformed 2nd-order cascaded transposed direct-form-II
%
% Usage: [y,V] = transpfilt(B,A,w0,x,V);
%            y = transpfilt(B,A,w0,x,V);         (default start with V=0)
%
% B,A = numerator and denominator Kx3 coefficient matrices
% w0  = center frequency in radians/sample
% x   = vector of input samples
% V   = Kx4 matrix of internal states
%
% y = vector of output samples
% V = updated internal states after the last x is processed
%
% notes: the row vector [V(i,1),V(i,2),V(i,3),V(i,4)] represents the internal states of the i-th
%        section, with V(i,1) being the closest to the output, and V(i,4), the farthest
%
%        implements Fig.5

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

function [y,V] = transpfilt(B,A,w0,x,V)

if nargin==0, help transpfilt; return; end

c0 = cos(w0); s0 = sqrt(1-c0^2);

if w0==0,    c0=1;  s0=0; end    % special cases
if w0==pi/2, c0=0;  s0=1; end
if w0==pi,   c0=-1; s0=0; end

K = size(B,1);

y = zeros(size(x));

if nargin==4, V = zeros(K,4); end

for n=1:length(x),
   x1 = x(n);
   for i=1:K,
      y1 = B(i,1)*x1 + V(i,1);
      w1 = B(i,2)*x1 - A(i,2)*y1 + V(i,3);
      V(i,1) = c0*w1 - s0*V(i,2); 
      V(i,2) = c0*V(i,2) + s0*w1;
      w2 = B(i,3)*x1 - A(i,3)*y1;
      V(i,3) = c0*w2 - s0*V(i,4);
      V(i,4) = c0*V(i,4) + s0*w2;
      x1 = y1;
   end
   y(n) = y1;
end



