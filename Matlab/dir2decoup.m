% dir2decoup.m - cascaded direct form to cascaded decoupled form
%
% Usage: [gamma,d] = dir2decoup(B,A)
%
% B,A = Kx3 matrices of numerator and denominator 2nd-order sections
%
% gamma = Kx2 matrix of reflection coefficients, columnwise, gamma = [gamma1,gamma2]
% d     = Kx3 matrix of ladder weights, columnwise, d = [d0,d1,d2]
%
% notes: 
%
%        it assumes that the first row B(1,:) is either 1st order or just a gain  
%
%          the transmission coefficients may be computed by tau = sqrt(1-gamma.^2) 

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

function [gamma,d] = dir2decoup(B,A)

if nargin==0, help dir2decoup; return; end

gamma1 = A(:,2)./(1+A(:,3));  tau1 = sqrt(1 - gamma1.^2);
gamma2 = A(:,3);              tau2 = sqrt(1 - gamma2.^2);

gamma = [gamma1,gamma2];

d = [1,A(1,2); A(1,2),1]\[B(1,1); B(1,2)];  d = [d', 0];    % first row

for i=2:size(B,1),
   d(i,:) = ([1, A(i,3), 1; A(i,2), A(i,2), 2; A(i,3), 1, 1]\[B(i,1); B(i,2); B(i,3)])';
   d(i,3) = d(i,3)/(tau1(i)*tau2(i));
end



   






 


