% cas2dir.m - cascade to direct form coefficients
%
% Usage: b = cas2dir(B)
%
% B = Kxp matrix of section coeffients (K sections of length p)
%
% b = length-(Kxp) direct-form coefficients
%
% notes: b is the convolution of the rows of B

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

function b = cas2dir(B)

if nargin==0, help cas2dir; return; end

K = size(B,1);

b = 1;

for i=1:K,
   b = conv(b,B(i,:));
end


