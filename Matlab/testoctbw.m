% testoctbw.m - test of the iterative solution of the bandwidth equation using octbw.m
%
% and plots the attenuation rate, Eq.(75), of the bandwidth iteration
% as a function of B and w0

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

iter = 8;

w0 = 0.4 * pi;  bw = 4;                    

for i=0:iter,                                 % calculate effective bw at each iteration
  [Dw,bweff(i+1),B] = octbw(w0,bw,i);
  err(i+1) = (bweff(i+1)-bw);                 % bw error
end

figure;

n=(0:iter);

fprintf('\n  n       bw       bw_error  \n');           % print results at each iteration
fprintf('----------------------------\n');
fprintf('%3.0f     %2.4f     %2.6f\n', [n', bweff', abs(err)']');

plot(n,log10(abs(err)),'r-', n,log10(abs(err)),'k.');   % note the exponential decrease

ylabel('log10(error)'); xlabel('iterations');
ylim([-15,1]); ytick(-14:2:0);
xtick(0:1:8); grid;

R = 2^(B/2); W0 = tan(w0/2);                           % use last B after iter iterations

a = (R^2+W0^2)*atan(W0/R)/((1+R^2*W0^2)*atan(W0*R))    % theoretical attenuation rate 

figure;
 
plot(n, err./(-a).^n, 'r-');           % test exponential decrease, this curve should be flat
                                       % theoretically, err(n) = const * (-a)^n
ylabel('error(n)/(-a)^n'); xlabel('iterations, n');
ylim([-2,2]); ytick(-2:1:2);
xtick(0:1:8); grid;


if 1

% ---------------------------------

N=5; type=0;   

G0 = 0; G = 12; GB = G-3;  GBoct = [GB,G,GB];

[w1,w2] = bandedge(w0,Dw);      % use the last Dw, after iter iterations 

w102 = log2([w1,w0,w2]/pi);

[B,A,Bh,Ah] = hpeq(N,G0,G,GB,w0,Dw,type);

a = 5; 
woct = linspace(-a, 0,1025); w = pi * 2.^woct;  

H = 20*log10(abs(fresp(B,A,w)));

figure;

plot(woct,H,'r-', w102, GBoct, 'b.');             % note the accuracy of the net octave bandwidth
                                                  % and the highly asymmetric band

ylabel('dB'); xlabel('octaves from Nyquist');

xlim([-a,0]); xtick(-a:1:0);
ylim([-1,15]); ytick(0:3:15);
grid;

end

% -------------------- surface plot of attenuation rate -----------
% -------------------- as a function of B and w0 -----------------


figure;

clear all;

[B,f0] = meshgrid(0:0.1:5, 0:1/40:1);

W = tan(pi*f0/2); 
R = 2.^(B/2);

a = (R.^2 + W.^2) .* atan(W./R) ./ (1+R.^2 .* W.^2) ./ atan(W.*R);

surf(B,f0,a); colorbar;

view(120,55);

xlabel('\it B'); ylabel('\it \omega_0/\pi'); 


