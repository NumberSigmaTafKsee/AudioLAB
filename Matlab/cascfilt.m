% cascfilt.m - filtering in cascade form (uses the built-in function filter)
%
% Usage: y = cascfilt(B,A,x)
%
% B,A = matrices whose rows are the numerator and denominator section coefficients
% x = vector of input samples
%
% y = vector of output samples

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
function y = cascfilt(B,A,x)

if nargin==0, help cascfilt; return; end

K = size(B,1);         % number of sections

y = x;

for i=1:K, 
   y = filter(B(i,:),A(i,:),y);
end


