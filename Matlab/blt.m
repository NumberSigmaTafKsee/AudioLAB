% blt.m - bilinear transformation of analog second-order sections
%
% Usage: [B,A,Bhat,Ahat] = blt(Ba,Aa,w0);
%
% Ba,Aa = Kx3 matrices of analog numerator and denominator coefficients (K sections)
% w0    = center frequency in radians/sample 
%
% B,A = Kx5 matrices of numerator and denominator coefficients (4th-order sections in z)
% Bhat,Ahat = Kx3 matrices of 2nd-order sections in the variable zhat
%
% notes: It implements the two-stage bilinear transformation: 
%                       s    -->    zhat    -->    z
%                  LP_analog --> LP_digital --> BP_digital
%
%        s = (zhat-1)/(zhat+1) = (z^2 - 2*c0*z + 1)/(z^2 - 1), with zhat = z*(c0-z)/(1-c0*z)
%
%        c0 = cos(w0), where w0 = 2*pi*f0/fs = center frequency in radians/sample
%
%        (B0 + B1*s + B2*s^2)/(A0 + A1*s + A2*s^2) = 
%        (b0h + b1h*zhat^-1 + b2h*zhat^-2)/(1 + a1h*zhat^-1 + a2h*zhat^-2) =
%        (b0 + b1*z^-1 + b2*z^-2 + b3*z^-3 + b4*z^-4)/(1 + a1*z^-1 + a2*z^-2 + a3*z^-3 + a4*z^-4)
%
%        column-wise, the input and output matrices have the forms:
%        Ba = [B0,B1,B2], Bhat = [b0h, b1h, b2h], B = [b0,b1,b2,b3,b4]
%        Aa = [A0,A1,A2], Ahat = [1,   a1h, a2h], A = [1, a1,a2,a3,a4]
        
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

function [B,A,Bhat,Ahat] = blt(Ba,Aa,w0)
    
if nargin==0, help blt; return; end

K = size(Ba,1);         % number of sections

B = zeros(K,5); A = zeros(K,5);
Bhat = zeros(K,3); Ahat = zeros(K,3); 

B0 = Ba(:,1); B1 = Ba(:,2); B2 = Ba(:,3);       % simplify notation
A0 = Aa(:,1); A1 = Aa(:,2); A2 = Aa(:,3);       % A0 may not be zero

c0 = cos(w0);       

if w0==0,    c0=1;  end;                        % make sure special cases are computed exactly
if w0==pi;   c0=-1; end;
if w0==pi/2; c0=0;  end;

i = find((B1==0 & A1==0) & (B2==0 & A2==0));    % find 0th-order sections (i.e., gain sections)

Bhat(i,1) = B0(i)./A0(i);
Ahat(i,1) = 1;

B(i,1) = Bhat(i,1);
A(i,1) = 1;

i = find((B1~=0 | A1~=0) & (B2==0 & A2==0));    % find 1st-order analog sections

D = A0(i)+A1(i);
Bhat(i,1) = (B0(i)+B1(i))./D;
Bhat(i,2) = (B0(i)-B1(i))./D;
Ahat(i,1) = 1;
Ahat(i,2) = (A0(i)-A1(i))./D;

B(i,1) = Bhat(i,1); 
B(i,2) = c0*(Bhat(i,2)-Bhat(i,1));
B(i,3) = -Bhat(i,2);
A(i,1) = 1;
A(i,2) = c0*(Ahat(i,2)-1);
A(i,3) = -Ahat(i,2);

i = find(B2~=0 | A2~=0);                        % find 2nd-order analog sections

D = A0(i)+A1(i)+A2(i);
Bhat(i,1) = (B0(i)+B1(i)+B2(i))./D;
Bhat(i,2) = 2*(B0(i)-B2(i))./D;
Bhat(i,3) = (B0(i)-B1(i)+B2(i))./D;
Ahat(i,1) = 1;
Ahat(i,2) = 2*(A0(i)-A2(i))./D;
Ahat(i,3) = (A0(i)-A1(i)+A2(i))./D;

B(i,1) = Bhat(i,1);
B(i,2) = c0*(Bhat(i,2)-2*Bhat(i,1));
B(i,3) = (Bhat(i,1)-Bhat(i,2)+Bhat(i,3))*c0^2 - Bhat(i,2);
B(i,4) = c0*(Bhat(i,2)-2*Bhat(i,3));
B(i,5) = Bhat(i,3);

A(i,1) = 1;
A(i,2) = c0*(Ahat(i,2)-2);
A(i,3) = (1-Ahat(i,2)+Ahat(i,3))*c0^2 - Ahat(i,2);
A(i,4) = c0*(Ahat(i,2)-2*Ahat(i,3));
A(i,5) = Ahat(i,3);

if c0==1 | c0==-1 	        % LP or HP shelving filter
   B = Bhat;                    % B,A are second-order
   A = Ahat;
   B(:,2) = c0*B(:,2);	        % change sign if w0=pi
   A(:,2) = c0*A(:,2);
   B(:,4:5) = 0;                % make them (K+1)x5
   A(:,4:5) = 0;                % for convenience in using fresp
end    







    
    
        











