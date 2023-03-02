% elliprf.m - elliptic rational function
%
% Usage: F = elliprf(N,k,w,tol)
%        F = elliprf(N,k,w)            (equivalent to tol=eps)
%
% N = analog filter order 
% k = ellptic modulus (defines transition band)
% w = vector of normalized frequencies
% tol = computational tolerance 
%
% F = function values at w
%
% Notes: F(w) = w^r * Prod((w^2 - z(i)^2)./(1 - w^2*k^2*z(i)^2) * (1 - k^2*z(i)^2)./(1 - z(i)^2))
%                i = 1:L, z(i) = cd(u(i)*K,k), u(i) = (2*i-1)/N
%
%        F(w) = cd(N*K1*u,k1), w = cd(u*K,k)
%
%        k1 = ellipdeg1(N,k) = solution of the degree equation, N*K'/K = K1'/K1

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

function F = elliprf(N,k,w,tol)

if nargin==0, help elliprf; return; end
if nargin==3, tol=eps; end

r = mod(N,2); L = (N-r)/2; 

i=1:L;  u = (2*i-1)/N; z = cde(u,k,tol);

F = ones(size(w));

if r==1,
   F = w;
end

for i=1:L,
   F = F .* (w.^2 - z(i)^2)./(1 - w.^2*k^2*z(i)^2) * (1 - k^2*z(i)^2)./(1 - z(i)^2);
end














