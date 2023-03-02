% ripple.m - calculate ripple factors from gains
%
% Usage: e = ripple(G,GB,G0)
%
% G  = peak/cut gain (all gains must be in dB, enter -Inf to get 0 in absolute units)
% GB = bandwidth gain
% G0 = reference gain 
%
% e = ripple factor, e = sqrt((G^2-GB^2)/(GB^2-G0^2))
%
% notes: must have G>GB>G0, or, G<GB<G0
%
%        G,GB,G0 are in dB and converted internally to absolute units
%
%        it carries out the following computations:
%
%        G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20);
%        e = sqrt((G^2-GB^2)/(GB^2-G0^2));

% --------------------------------------------------------------
% Copyright (c) 2005 by Sophocles J. Orfanidis
% 
% Address: Sophocles J. Orfanidis                       
%          ECE Department, Rutgers University          
%          94 Brett Road, Piscataway, NJ 08854-8058
% Email:   orfanidi@ece.rutgers.edu
% Date:    June 15, 2005
% 
% Reference: Sophocles J. Orfanidis, "High-Order Digital Parametric Equalizer 
%            Design," J. Audio Eng. Soc., vol.53, pp. 1026-1046, November 2005.
%
% Web Page: http://www.ece.rutgers.edu/~orfanidi/hpeq
% 
% tested with MATLAB R11.1 and R14
% --------------------------------------------------------------

function e = ripple(G,GB,G0)

if nargin==0, help ripple; return; end

G0 = 10^(G0/20); G = 10^(G/20); GB = 10^(GB/20);

e = sqrt((G^2-GB^2)/(GB^2-G0^2));    

