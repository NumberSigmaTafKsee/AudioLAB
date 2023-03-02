% fresp_a.m - frequency response of cascaded analog filter 
%
% Usage: H = fresp_a(B,A,f,f0)
%        H = fresp_a(B,A,f)       (equivalent to f0=0)
%
% B = K x p matrix of K length-p numerator sections
% A = L x q matrix of L length-q denominator sections 
% f = vector of frequencies in Hz, (row or column)
% f0 = center frequency in Hz, (default f0=0)
%
% H = vector of frequency responses, has same size as f
%
% notes: B,A are extended internally to sizes Mxp and Mxq, respectively, where M = max(K,L)
%
%        B,A of a single section (K=L=1) must be entered as rows
%
%        it evaluates the product H(w) = [B1(w)/A1(w)]*[B2(w)/A2(w)]*...*[BM(w)]/[AM(w)]
%
%        optionally, it applies the lowpass-to-bandpass analog transformation
%        s => s + w0^2/s, to recenter a lowpass frequency response at f=f0
%
%        it uses the built-in function POLYVAL

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

function H = fresp_a(B,A,f,f0)

if nargin==0, help fresp_a; return; end
if nargin==3, f0=0; end

w0 = 2*pi*f0;

K = size(B,1);		           % number of numerator sections
L = size(A,1);		           % number of denominator sections
M = max(K,L);                      % extend to a common number of sections

B(K+1:M,:)=0; B(K+1:M,1)=1;        % extra sections are trivial, i.e., [1,0,0,...,0]
A(L+1:M,:)=0; A(L+1:M,1)=1;

s = 2*pi*j*f;

if f0==0,
   s1 = s;
elseif f0==Inf;
   s1 = 1./s;
else
   s1 = s + w0^2./s;
end

H = 1;

for i=1:M,
   H = H .* polyval(fliplr(B(i,:)), s1) ./ polyval(fliplr(A(i,:)), s1);	
end




