% ellipdeg.m - solves the degree equation in analog elliptic filter design
%
% Usage: k = ellipdeg(N,k1,tol);  
%        k = ellipdeg(N,k1,M);    (M = fixed number of Landen iterations)
%        k = ellipdeg(N,k1)       (uses tol=eps)
%
% N = analog filter order
% k1 = elliptic modulus for stopband band, that is, k1 = ep/ep1
% tol = tolerance, e.g., tol=1e-10, default is tol = eps
%
% M = uses a fixed number of Landen iterations, typically, M = 4-5
%
% k = elliptic modulus for transition band, that is, k = WB/W1
%
% Notes: solves the degree equation N*K'/K = K1'/K1 for k in terms of N,k1
%
%        it uses Jacobi's exact solution k' = (k1')^N * (Prod_i(sn(ui*K1',k1')))^4
% 
%        when k1 is very small, it uses the function ELLIPDEG2 to avoid numerical
%        problems arising from the computation of sqrt(1-k1^2)
%
%        to solve for k1, given N,k, use the function ELLIPDEG1

% -------------------------------------------------------------------------
% Copyright (c) 2005 by Sophocles J. Orfanidis
% 
% Address: Sophocles J. Orfanidis                       
%          ECE Department, Rutgers University          
%          94 Brett Road, Piscataway, NJ 08854-8058, USA
%
% Email:   orfanidi@ece.rutgers.edu
% Date:    June 15, 2005, revised February 19, 2006
% 
% Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
%
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% 
% tested with MATLAB R11.1 and R14
% -------------------------------------------------------------------------

function k = ellipdeg(N,k1,tol)

if nargin==0, help ellipdeg; return; end
if nargin==2, tol=eps; end

L = floor(N/2); 

i = 1:L;  ui = (2*i-1)/N;

kmin = 1e-6;

if k1 < kmin
   k = ellipdeg2(1/N,k1,tol);              % use nome method when k1 is too small
else                                       
   kc = sqrt(1-k1^2);			   % complement of k1 
   kp = kc^N * (prod(sne(ui,kc)))^4;       % complement of k
   k = sqrt(1-kp^2); 
end


% when k1 is near 1, then so is k, because k1 < k < 1, and if the calculated kp 
% becomes less than about eps, then k = sqrt(1-kp^2) would be inaccurate. 
% To avoid this, N may not be too large. The maximum usable value of N consistent
% with Matlab's floating point accuracy is when k becomes equal to about k=1-eps,
% which gives Nmax = log(q1)/log(q), with q1 = nome(k1) and q = nome(1-eps), or,
% approximately, Nmax = -3.86 * log(q1) 









