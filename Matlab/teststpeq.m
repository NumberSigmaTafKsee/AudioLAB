% teststpeq.m - testing of the biquad state-space form

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

clear all;

fs=40; 
f0 = 14;  Df = 2;  w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;

c0 = cos(w0); s0 = sqrt(1-c0^2);

if w0==0,    c0==1; s0=0; end   % special cases
if w0==pi/2, c0=0;  s0=1; end
if w0==pi,   c0=-1; s0=0; end

sig = (c0>=0)-(c0<0);           % sign of c0, but sign(0) = 1

G0db = 0; Gdb = 12; GBdb = 9;

G0 = 1; G = 10^(Gdb/20); GB = 10^(GBdb/20); 

e = sqrt((G^2-GB^2)/(GB^2-G0^2)); be = tan(Dw/2)/e;


b1 = [G0+G*be, -2*G0*c0, G0-G*be]/(1+be);
a1 = [1, -2*c0/(1+be), (1-be)/(1+be)];

[b2,a2] = hpeq(1,G0db,Gdb,GBdb,w0,Dw,0);

A = [c0, be+s0; be-s0, c0]/(1+be); 
B = sqrt(2*be) * [sqrt(1-s0); -sig*sqrt(1+s0)]/(1+be);
C = 0.5*sqrt(2*be)*(G-G0)*[sig*sqrt(1+s0), -sqrt(1-s0)]/(1+be);
D = (G0+G*be)/(1+be);

[A1,B1,C1,D1] = stpeq(G0db, Gdb, GBdb, w0, Dw);

ABCD_diff = norm(A-A1)+norm(B-B1)+norm(C-C1)+norm(D-D1)

figure;

f = linspace(0,20,501); w = 2*pi*f/fs; z = exp(j*w);

for i=1:length(z),
  H(i) = D + C*((z(i)*eye(2)-A)\B);   % H(z) = D + C * inv((z*I - A)) * B
end

H2 = fresp(b2,a2,w);

H_diff = norm(H-H2)

plot(f,20*log10(abs(H)),'b-', f,20*log10(abs(H2)),'--r');

ylim([-14 14]); ytick(-12:3:12);
xlim([0,20]); xtick(0:2:20);
grid;


