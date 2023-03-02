% srem.m - symmetrized rem
%
% Usage: Z = srem(X,Y)
%
% X = real-valued vector
% Y = positive scalar
%
% Z = has same size as X, and lies in the interval [-Y/2, Y/2]
%
% Notes: same syntax as REM, but it brings the result Z into the 
%        symmetric interval [-Y/2, Y/2]
%
%        Z = rem(X,Y)
%        Z = Z - Y.*sign(Z).*(abs(Z)>Y/2)

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

function Z = srem(X,Y)

if nargin==0, help srem; return; end

Z = rem(X,Y);                       % bring into interval [-Y,Y]

Z = Z - Y.*sign(Z).*(abs(Z)>Y/2);   % bring into interval [-Y/2,Y/2]


