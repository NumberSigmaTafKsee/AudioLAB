% dne.m - dn elliptic function with normalized real argument
%
% Usage: w = dne(u,k,tol)  (e.g., tol=1e-8)
%        w = dne(u,k,M)    (M=integer)
%        w = dne(u,k)      (equivalent to tol=eps)
%
% u = arbitrary vector of real numbers
% k = elliptic modulus (0 <= k < 1)
% tol = tolerance, e.g., tol=1e-8, default is tol = eps
%
% M = use a fixed number of Landen iterations, typically, M = 4-5
%
% w = the value of dn(u*K,k), w has the same size as u
%
% Notes: u is in units of the quarterperiod K, thus, dn(x,k) = dne(x/K,k)
%
%        K = K(k), K' = K'(k) = K(k'), k' = sqrt(1-k^2)
%
%        it uses SNE and the formula dn(x,k) = sqrt(1 - k^2 * sn(x,k)^2)
        
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

function w = dne(u,k,tol)

if nargin==0, help dne; return; end
if k==1, disp('k may not be equal to 1'); return; end
if nargin==2, tol=eps; end
if ~isreal(u), fprintf('\n u must be real \n\n'); return; end

w = sqrt(1 - k^2 * sne(u,k,tol).^2);





