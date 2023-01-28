% ellipdeg1.m - solves the degree equation in analog elliptic filter design
%
% Usage: k1 = ellipdeg1(N,k,tol);  
%        k1 = ellipdeg1(N,k,M);    (M = fixed number of Landen iterations)
%        k1 = ellipdeg1(N,k)       (uses tol=eps)
%
% N = analog filter order
% k = elliptic modulus for transition band, that is, k = WB/Ws
% tol = tolerance, e.g., tol=1e-10, default is tol = eps
%
% M = uses a fixed number of Landen iterations, typically, M = 4-5
%
% k1 = elliptic modulus for stopband, k1 = e/es
%
% Notes: solves the degree equation N*K'/K = K1'/K1 for k1 in terms of N,k
%

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

function k1 = ellipdeg1(N,k,tol)

if nargin==0, help ellipdeg1; return; end
if nargin==2, tol=eps; end

L = floor(N/2); 

i = 1:L;  ui = (2*i-1)/N;

k1 = k^N * (prod(sne(ui,k)))^4;       










