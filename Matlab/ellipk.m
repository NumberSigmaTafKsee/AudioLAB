% ellipk - complete elliptic integral of first kind
% 
% Usage: [K,Kprime] = ellipk(k,tol)
%        [K,Kprime] = ellipk(k)          (equivalent to tol=eps)
%        [K,Kprime] = ellipk(k,M)        (equivalent to using M Landen iterations)
%
% k   = elliptic modulus 
% tol = tolerance, e.g., tol = 10^(-10), default value is tol = eps = 2.22e-16 
% tol = M = fixed number of Landen iterations resulting in length(v) = M
%
% K      = quarter period K(k)
% Kprime = quarter period K'(k) = K(k'), k' = sqrt(1-k^2)
% 
% Notes: first it constructs the Landen vector of descending moduli, v = landen(k), 
%        and then it computes K = prod(1+v)) * pi/2
%
%        when k is near 0, it uses the approximation Kprime = -log(k/4)
%        when k is near 1, it uses the approximation K = -log(k'/4)
%        however "near 1" means that k cannot be nearer to 1 than eps
%
%        the nome is q = exp(-pi*K'/K) 
%
%        produces the same answer as the built-in function ELLIPKE, K = ellipke(k^2)

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

function [K,Kprime] = ellipk(k,tol)

if nargin==0, help ellipk; return; end
if nargin==1, tol=eps; end

kmin = 1e-6; 
kmax = sqrt(1-kmin^2);

if k==1
   K = Inf;
elseif k > kmax                   % floating point resolution limits this to kmax < k < 1-eps
   kp = sqrt(1-k^2);              % k > kmax is equivalent to kp < kmin
   L = -log(kp/4);
   K = L + (L-1)*kp^2/4;                          
   % K = L + (L-1)*kp^2/4 + (L-7/6)*9*kp^4/64;    
else
   v = landen(k,tol);  
   K = prod(1+v) * pi/2;
end

if k==0
   Kprime = Inf;
elseif k < kmin                                      
   L = -log(k/4);
   Kprime = L + (L-1)*k^2/4;                           % O(k^4) terms do not affect floating point accuracy 
   % Kprime = L + (L-1)*k^2/4 + (L-7/6)*9*k^4/64;      % O(k^4) terms would be needed if kmin = (eps)^(1/4)                  
else                                   
   kp = sqrt(1-k^2);                             
   vp = landen(kp,tol);  
   Kprime = prod(1+vp) * pi/2;
end





