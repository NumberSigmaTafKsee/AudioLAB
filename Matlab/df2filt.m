% df2filt.m - filtering in frequency-shifted 2nd-order cascaded direct form II
%
% Usage: [y,U] = df2filt(B,A,w0,x,U);
%            y = df2filt(B,A,w0,x);          (default start with U=0)
%
% B,A = numerator and denominator Kx3 coefficient matrices
% w0  = recentered frequency in radians/sample
% x   = vector of input samples
% U   = Kx4 matrix of internal states
%
% y = vector of output samples
% U = updated internal states after the last x is processed
%
% notes: the row vector [U(i,1),U(i,2),U(i,3),U(i,4)] represents the internal states of the i-th
%        section, with U(i,1) being the closest to the output, and U(i,4), the farthest
%
%        implements the trasposed of Fig.5

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

function [y,U] = df2filt(B,A,w0,x,U)

if nargin==0, help df2filt; return; end

c0 = cos(w0); s0 = sqrt(1-c0^2);

if w0==0,    c0=1;  s0=0; end    % special cases
if w0==pi/2, c0=0;  s0=1; end
if w0==pi,   c0=-1; s0=0; end

K = size(B,1);

y = zeros(size(x));     % makes y a column or row according to x

if nargin==4, U = zeros(K,4); end

for n=1:length(x),
   y1 = x(n);
   for i=1:K,
      v1 = c0*U(i,1) + s0*U(i,2);
      v2 = c0*U(i,3) + s0*U(i,4);
      v0 = y1 - A(i,2)*v1 - A(i,3)*v2;
      y1 = B(i,1)*v0 + B(i,2)*v1 + B(i,3)*v2;    
      U(i,4) = c0*U(i,4) - s0*U(i,3);
      U(i,3) = v1;
      U(i,2) = c0*U(i,2) - s0*U(i,1);
      U(i,1) = v0;
   end
   y(n) = y1;
end



