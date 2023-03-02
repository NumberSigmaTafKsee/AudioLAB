% stpeq.m - minimum-noise state-space realization of biquadratic digital parametric equalizer
%
% Usage: [A,B,C,D] = stpeq(G0, G, GB, w0, Dw); 
%
% G0 = reference gain in dB (all gains must be in dB, enter G0=-inf to get G0=0 in absolute units)
% G  = peak/cut gain in dB
% GB = bandwidth gain in dB
% w0 = peak/cut center frequency in units of [radians/sample], i.e., w0=2*pi*f0/fs
% Dw = bandwidth in [radians/sample]
%
% A,B,C,D = two-dimensional ABCD state-space parameters (one-dimenslional if w0=0 or w0=pi)
%
% notes: w0 = 2*pi*f0/fs, Dw = 2*pi*Df/fs, with f0,Df,fs in Hz

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

function [A,B,C,D] = stpeq(G0, G, GB, w0, Dw)

if nargin==0, help stpeq; return; end

G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20);   % gains in absolute units

e = sqrt((G^2-GB^2)/(GB^2-G0^2));
b = tan(Dw/2)/e;

c0 = cos(w0); s0 = sqrt(1-c0^2);

if w0==0,    c0=1;  s0=0; end        % special cases
if w0==pi,   c0=-1; s0=0; end
if w0==pi/2, c0=0;  s0=1; end

p = (c0>=0) - (c0<0);                % sign of c0, but p=1 if c0=0

A = [c0, b+s0; b-s0, c0]/(1+b); 
B = sqrt(2*b)*[sqrt(1-s0); -p*sqrt(1+s0)]/(1+b);
C = sqrt(2*b)*(G-G0)*[p*sqrt(1+s0), -sqrt(1-s0)]/(1+b)/2; 
D = (G0+G*b)/(1+b);

if c0==1 | c0==-1,                            % first-order shelving filters
   A = [c0*(1-b)/(1+b), 0; 0, 0];
   B = 2*sqrt(b)/(1+b) * [1; 0];
   C = c0*sqrt(b)*(G-G0)/(1+b) * [1,0];
   D = (G0+G*b)/(1+b);
end




