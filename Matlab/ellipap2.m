% ellipap2.m - analog lowpass elliptic filter design
%
% Usage: [z,p,H0,B,A] = ellipap2(N,Ap,As)
%
% N  = filter order
% Ap = passband attenuation in dB
% As = stopband attenuation in dB
%
% z  = vector of normalized filter zeros (in units of the passband frequency Wp = 2*pi*fp)
% p  = vector of normalized filter poles
% H0 = DC gain factor
% B  = matrix whose rows are the first- and second-order numerator coefficients
% A  = matrix whose rows are the first- and second-order denominator coefficients
%
% notes: serves as a substitute for the built-in function ELLIPAP
%
%        the gain factor g returned by ELLIPAP is related to the dc gain by g = abs(H0*prod(p)/prod(z))
%
%        N = 2*L+r, r = mod(N,2), L = floor(N/2) = no. second-order sections
%
%        length(p) = N, length(z) = 2*L
%
%        transfer function: H(s) = H0 * [1/(1-s/p0)]^r * Prod_{i=1}^L [(1-s/zi)(1-s/zi^*)]/[(1-s/pi)(1-s/pi^*)]
%
%        normalized s-plane variable: s = j*Om/Om_p

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

function [z,p,H0,B,A] = ellipap2(N,Ap,As)

if nargin==0, help ellipap2; return; end

Gp = 10^(-Ap/20);                                      % passband gain
ep = sqrt(10^(Ap/10) - 1);                             % ripple factors
es = sqrt(10^(As/10) - 1);

k1 = ep/es;                                 
k = ellipdeg(N,k1);                                    % solve degree equation

L = floor(N/2); r = mod(N,2);                          % L is the number of second-order sections
i = (1:L)'; ui = (2*i-1)/N; zeta_i = cde(ui,k);        % zeros of elliptic rational function

z = j./(k*zeta_i);                                     % filter zeros = poles of elliptic rational function

v0 = -j*asne(j/ep, k1)/N;                              % solution of sn(j*v0*N*K1,k1) = j/ep

p = j*cde(ui-j*v0, k);                                 % filter poles
p0 = j*sne(j*v0, k);                                   % first-order pole, needed when N is odd

B = [ones(L,1), -2*real(1./z), abs(1./z).^2];          % second-order numerator sections
A = [ones(L,1), -2*real(1./p), abs(1./p).^2];          % second-order denominator sections

if r==0,                                               % prepend first-order sections
  B = [Gp, 0, 0; B];                                   % dc gain is Gp for N even
  A = [1, 0, 0; A];
else
  B = [1, 0, 0; B];                                    % dc gain is 1 for N odd
  A = [1, -real(1/p0), 0; A];
end

z = cplxpair([z; conj(z)]);                            % append conjugate zeros
p = cplxpair([p; conj(p)]);                            % append conjugate poles
if r==1, p = [p; p0]; end                              % append first-order pole when N is odd

H0 = Gp^(1-r);                                         % dc gain




