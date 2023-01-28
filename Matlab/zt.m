% zt.m - evaluates z-transform of a finite-duration sequence
%
% Usage: A = zt(a, z);
%
% a = order-M filter, a = [a0,a1,...,aM], (entered as column or row)
% z = any vector of non-zero complex numbers, (entered as column or row)
%
% A = polynomial A(z) evaluated at the matrix of z's, same size as z
%
% notes: uses Hoerner's rule
%
%        has the same output as A = polyval(flip(a),1./z)

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

function A = zt(a, z)

if nargin==0, help zt; return; end
if ~isempty(find(z==0)), disp('z must not have any zero entries'); return; end

[K,L] = size(z);

a = a(:);                     %  make a,z into columns
z = z(:);                     

z = 1./z;                           

A = 0;                          
for n = length(a):-1:1,
    A = a(n) + z .* A;       
end

A = reshape(A,K,L);

