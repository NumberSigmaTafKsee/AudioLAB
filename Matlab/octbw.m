% octbw.m - octave to linear bandwidth computation
%
% Usage: [Dw,bn,B] = octbw(w0,b,iter);
%        [Dw,bn,B] = octbw(w0,b);        (equivalent to iter=0)
%
% w0   = center frequency in rads/sample (may not be 0 or pi)
% b    = desired bandwidth in octaves, i.e., b = log2(w2/w1)
% iter = number of iterations to solve the octave bandwidth equation (typically, 0-3)
%
% Dw = linear bandwidth, Dw = w2 - w1
% bn = calculated bandwidth in octaves (after n=iter iterations)
% B  = calculated prewarped analog bandwidth
%
% notes: Starting with Bristow-Johnson's approximation B = w0/sin(w0) * b,
%        it iteratively solves (for B in terms of b) the octave bandwidth equation:
%
%        atan(2^(B/2)*t0) / atan(2^(-B/2)*t0) = 2^b, where t0=tan(w0/2)
%
%        After n=iter of iterations, which update the value of B, the octave bandwidth bn is 
%        computed from:  2^bn = atan(2^(B/2)*t0)/atan(2^(-B/2)*t0), 
%        and the corresponding linear bandwidth from: tan(Dw/4) = sin(w0) * sinh(B*log(2)/2)
%
%        iter=0 corresponds to the Bristow-Johnson approximation
%
%        typically, 1-3 iterations are sufficient for all w0 across 0 < w0 < 0.9*pi and 
%        any reasonable, even very large, value of b, with the approximation error 
%        abs(bn-b) essentially decreasing exponentially with the number of iterations
%
%        once the quantities w0,Dw are known, the bandedge frequencies are computed by
%        [w1,w2] = bandedge(w0,Dw), and the EQ filter design can proceed using HPEQ
%
%        the computed bandedge frequencies w1,w2 always lie within the Nyquist interval (0,pi),
%        and although they are b-octaves apart (actually, bn octaves), they do not lie 
%        symmetrically at +b/2 and -b/2 about w0, and may result in a very asymmetric band, 
%        especially for large w0's  
%
% Reference: R. Bristow-Johnson, "The Equivalence of Various Methods of Computing Biquad 
%            Coefficients for Audio Parametric Equalizers," 97th Convention of the AES, 
%            San Francisco, November 1994, {\em AES Preprint 3906}

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

function [Dw,bn,B] = octbw(w0,b,iter)

if nargin==0, help octbw; return; end
if nargin==2, iter=0; end

B = w0/sin(w0) * b;                         % Bristow-Johnson's approximation

t0 = tan(w0/2);                             

r = 2^(b);                                  % desired ratio of bandedge frequencies, r = w2/w1
R = 2^(B/2);                                % introduced for notational convenience

for n=1:iter,                               % iterative solution of atan(t0/R) = atan(R*t0) / r
   R = t0 / tan(atan(R*t0)/r);        
end

B = 2*log2(R);                              % calculated analog bandwidth

bn = log2(atan(R*t0)/atan(t0/R));           % calculated bandwidth in octaves, error = abs(bn-b)

Dw = 2*atan(sin(w0) * sinh(B*log(2)/2));    % calculated linear bandwidth in rads/sample




   





