% acde.m - inverse of cd elliptic function
%
% Usage: u = acde(w,k,tol)  (e.g., tol=1e-8)
%        u = acde(w,k,M)    (M=integer)
%        u = acde(w,k)      (equivalent to tol=eps)
%
% w = arbitrary vector of complex numbers on the u-plane
% k = elliptic modulus (0 <= k < 1)
% tol = tolerance, e.g., tol=1e-8, default is tol = eps
% M = use a fixed number of Landen iterations, typically, M = 4-5
%
% u = the solution of cd(u*K,k) = w, u has the same size as w, u=1 if w=0
%
% Notes: u is in units of the quarterperiod K, so that w = cd(u*K,k) = cde(u,k) is inverted
%        by  u = acde(w,k)
%
%        K = K(k), K' = K'(k) = K(k'), k' = sqrt(1-k^2), R=K'/K                    
%                                                                                  4 | 3 
%        u is reduced into the rectangle, 0<Re(u)<2, -R<Im(u)<R, quadrant mapping  ------, centered at u=1
%                                                                                  1 | 2 
%        k may not be 1, because K=inf, but note that cd(u,1) = 1                  
%        
%        it uses the Landen/Gauss transformation of descending moduli
%        to reduce the problem to cos(u*pi/2) = w, with solution u = 2*acos(w)/pi
%
%        the tolerance is that of computing the Landen vector v = landen(k,tol)
%
%        CDE and ACDE are inverses of each other
%        see also LANDEN, SNE, ASNE, ELLIPK, ELLIPDEG, and the built-in functions ELLIPJ, ELLIPKE
%
%        it uses the symmetrized REM function, SREM
        
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

function u = acde(w,k,tol)

if nargin==0, help acde; return; end
if k==1, disp('k may not be equal to 1'); return; end
if nargin==2, tol=eps; end

v = landen(k,tol);

for n = 1:length(v),                   	              % descending Landen/Gauss transformations 
    if n==1, v1=k; else v1=v(n-1); end                % v1 stands for v(n-1) = 2*sqrt(v(n))/(1+v(n))
    w = w./(1 + sqrt(1 - w.^2 * v1^2)) * 2/(1+v(n));  % solves 1/w(n-1) = (1/w(n) + v(n)*w(n))/(1+v(n)), for 1/w(n)
end                                                   % where w(n) = cd(u*K(n),v(n)) = cde(u,v(n));

u = 2/pi * acos(w);                                   % gives positive Re(u)

u(find(w==1))=0;				      % set exactly u=0 when w=1

[K,Kprime] = ellipk(k,tol);  R = Kprime/K; 

u = srem(real(u),4) + j*srem(imag(u),2*R);            % -R<Im(u)<R, and because Re(u)>0 ==> 0<Re(u)<2

