% dir2state.m - direct to optimum state-space form for second-order filters
%
% Usage: [A,B,C,D] = dir2state(b,a)
%
% b,a = Kx3 matrices of numerator and denominator 2nd-order sections
%
% A = 2x2xK array of state-transition matrices
% B = 2xK array 
% C = Kx2 array
% D = 1xK array    
% 
% notes: the ABCD parameters for the i-th section are A(:,:,i), B(:,i), C(i,:), D(i)
%
%        used in conjunction with STATEFILT for filtering in state-space form
%
%        Reference: C. W. Barnes, On the Design of Optimal State-Space Realizations of Second-Order 
%                   Digital Filters, IEEE Trans. Circuits Syst., CAS-31, 602 (1984).

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

function [A,B,C,D] = dir2state(b,a)

if nargin==0, help dir2state; return; end

K = size(b,1);    % number of second-order sections

if K==1,                     % ABCD parameters for single 2nd-order section

   D = b(1);

   q1 = b(2) - a(2)*b(1); q2 = b(3) - a(3)*b(1);

   if b(2)==0 & a(2)==0 & b(3)==0 & a(3)==0,           % 0th order - just a gain
      A = [0, 0; 0, 0]; 
      B = [0; 0]; 
      C = [0, 0];    
   elseif b(3)==0 & a(3)==0,                           % 1st order
      b1 = sqrt(1-a(2)^2); c1 = q1/b1;
      A = [-a(2), 0; 0, 0]; 
      B = [b1; 0]; 
      C = [c1, 0]; 
   else                                                % 2nd order
      s = -a(2)/2; w = sqrt(a(3) - s^2); p = s+j*w;
      ar = q1/2; ai = -(q1*s+q2)/(2*w); alpha = ar+j*ai;
      P = abs(alpha)/(1-abs(p)^2); Q = imag(alpha/(1-p^2)); k = sqrt((P+Q)/(P-Q));

      if ar~=0,
         b1 = sqrt((abs(alpha)-ai)/(P-Q)); b2 = -sqrt((abs(alpha)+ai)/(P+Q))*sign(ar);
         c1 = ar/b1; c2 = ar/b2;
      else
         if ai>0,
            b1 = 0; b2 = -sqrt(2*abs(alpha)/(P+Q));
            c1 = sqrt(2*abs(alpha)*(P-Q)); c2 = 0;
         else
            b1 = sqrt(2*abs(alpha)/(P-Q)); b2 = 0;
            c1 = 0; c2 = -sqrt(2*abs(alpha)*(P+Q));
         end
      end

      A = [s, w*k; -w/k, s];
      B = [b1; b2]; 
      C = [c1,c2];
   end

else                         % recursive definition

   for i=1:K, 
      [A(:,:,i),B(:,i),C(i,:),D(i)] = dir2state(b(i,:),a(i,:));
   end

end





