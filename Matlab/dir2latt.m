% dir2latt.m - cascaded direct form to cascaded normalized lattice form
%
% Usage: [gamma,d] = dir2latt(B,A)
%
% B,A = Kx3 matrices of numerator and denominator 2nd-order sections
%
% gamma = Kx2 matrix of reflection coefficients, columnwise, gamma = [gamma1,gamma2]
% d     = Kx3 matrix of ladder weights, columnwise, d = [d0,d1,d2]
%
% example: B = [10.0000         0         0      A = [1.0000         0         0
%                2.0000    3.0000         0           1.0000    0.8000         0
%                4.0000    5.0000    6.0000           1.0000   -0.3000   -0.4000
%                7.0000    8.0000    9.0000]          1.0000         0   -0.9000]
% 
%          gamma = [0.0   0.0   tau = [1.000  1.0000   d = [10.0000   0.0000   0
%                   0.8   0.0          0.600  1.0000        -0.6667   3.0000   0
%                  -0.5  -0.4          0.866  0.9165        12.3468   7.4194   6
%                   0.0  -0.9]         1.000  0.4359]       34.6418  18.3533   9]
%
%          gain factors and 1st-order sections are also viewed as 2nd-order sections
%
%          the transmission coefficients may be computed by tau = sqrt(1-gamma.^2) 
%
%          see Fig.6

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

function [gamma,d] = dir2latt(B,A)

if nargin==0, help dir2latt; return; end

gamma1 = A(:,2)./(1+A(:,3));  tau1 = sqrt(1 - gamma1.^2);
gamma2 = A(:,3);              tau2 = sqrt(1 - gamma2.^2);

d2 = B(:,3);
d1 = (B(:,2) - A(:,2).*d2)./tau2;
d0 = (B(:,1) - gamma1.*tau2.*d1 - A(:,3).*d2)./(tau1.*tau2);

gamma = [gamma1,gamma2];

d = [d0,d1,d2];




