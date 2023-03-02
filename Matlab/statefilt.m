% statefilt.m - filtering in frequency-shifted cascaded state-space form
%
% Usage: [y,S] = statefilt(A,B,C,D,w0,x,S);
%            y = statefilt(A,B,C,D,w0,x);        (default start with S=0)
%
% A  = 2x2xK array of state-transition matrices
% B  = 2xK array 
% C  = Kx2 array
% D  = 1xK array    
% w0 = center frequency (rads/sample)
% x  = vector of input samples
% S  = 2x2xK matrix of input internal states
% 
% y = vector of output samples
% S  = 2x2xK matrix of output internal states, after processing last x
%
% notes: for the i-th section the ABCD parameters are: A(:,:,i), B(:,i), C(i,:), D(i)
%        and the two 2-dimensional internal states:    S(:,1,i), S(:,2,i), 
%        with S(:,1) being clser to the output y
%
%        may be used in conjunction with DIR2STATE to calculate the ABCD parameters

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

function [y,S] = statefilt(A,B,C,D,w0,x,S)

if nargin==0, help statefilt; return; end

K = size(A,3);                     % number of sections

c0 = cos(w0); s0 = sqrt(1-c0^2);

if w0==0,    c0=1;  s0=0; end      % special cases
if w0==pi/2, c0=0;  s0=1; end
if w0==pi,   c0=-1; s0=0; end

if nargin==6, S = zeros(2,2,K); end

y = x;

for n=1:length(x),
  for i=1:K,
    U = A(:,:,i) * S(:,1,i) + B(:,i) * y(n);                  
    y(n) = C(i,:) * S(:,1,i) + D(i) * y(n);
    S(:,1,i) = c0 * U - s0 * S(:,2,i);
    S(:,2,i) = s0 * U + c0 * S(:,2,i);
  end
end





