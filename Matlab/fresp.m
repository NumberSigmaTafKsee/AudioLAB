% fresp.m - frequency response of a cascaded IIR filter at a frequency vector w
%
% Usage: H = fresp(B,A,w,w0)
%        H = fresp(B,A,w)       (equivalent to w0=0)
%
% B = K x p matrix of K length-p numerator sections
% A = L x q matrix of L length-q denominator sections 
% w = vector of frequencies in rads/sample, (row or column)
% w0 = center frequency in rads/sample, (default w0=0)
%
% H = vector of frequency responses, has same size as w
%
% notes: B,A are extended internally to sizes Mxp and Mxq, respectively, where M = max(K,L)
%
%        B,A of a single section (K=L=1) must be entered as rows
%
%        it evaluates the product H(w) = [B1(w)/A1(w)]*[B2(w)/A2(w)]*...*[BM(w)]/[AM(w)]
%
%        optionally, it applies the lowpass-to-bandpass bilinear transformation
%        zhat = z.*(c0-z)./(1-c0*z), c0=cos(w0), to recenter a lowpass frequency 
%        response at w=w0
%
%        it uses the z-transform evaluation function ZT

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


function H = fresp(B,A,w,w0)

if nargin==0, help fresp; return; end
if nargin==3, w0=0; end

K = size(B,1);		           % number of numerator sections
L = size(A,1);		           % number of denominator sections
M = max(K,L);                      % extend to a common number of sections

B(K+1:M,:)=0; B(K+1:M,1)=1;        % extra sections are trivial, i.e., [1,0,0,...,0]
A(L+1:M,:)=0; A(L+1:M,1)=1;

z = exp(j*w);

if w0==0,                          % special cases
   zhat = z;
elseif w0==pi/2,
   zhat = -z.^2;
elseif w0==pi,
   zhat = -z;
else
   c0 = cos(w0);
   zhat = z.*(c0-z)./(1-c0*z);     % lowpass-to-bandpass transformation
end

H = 1;

for i=1:M,
   H = H .* zt(B(i,:),zhat) ./ zt(A(i,:),zhat);	
end




