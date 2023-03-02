

function [y] = combFilter(x, a, b, z)
L = length(x);
y = zeros(L, 1);

for i = 1:z
    y(i) = x(i);
end

% filtering
for i = z+1:L
    y(i) = x(i) + a * x(max(1, i-z)) + b * y(max(1, i-z));  
end
*/

% Authors: P. Dutilleux, U Zölzer
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------


x=zeros(100,1);x(1)=1; % unit impulse signal of length 100
g=0.5;
Delayline=zeros(10,1);% memory allocation for length 10
for n=1:length(x);
	y(n)=x(n)+g*Delayline(10);
	Delayline=[x(n);Delayline(1:10-1)];
end;


% Authors: P. Dutilleux, U Zölzer
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

x=zeros(100,1);x(1)=1; % unit impulse signal of length 100
g=0.5;
Delayline=zeros(10,1); % memory allocation for length 10
for n=1:length(x);
	y(n)=x(n)+g*Delayline(10);
	Delayline=[y(n);Delayline(1:10-1)];
end;



function y = apbandpass (x, Wc, Wb)
% y = apbandpass (x, Wc, Wb)
% Author: M. Holters
% Applies a bandpass filter to the input signal x.
% Wc is the normalized center frequency 0<Wc<1, i.e. 2*fc/fS.
% Wb is the normalized bandwidth 0<Wb<1, i.e. 2*fb/fS.
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

c = (tan(pi*Wb/2)-1) / (tan(pi*Wb/2)+1);
d = -cos(pi*Wc);
xh = [0, 0];
for n=1:length(x)
  xh_new = x(n) - d*(1-c)*xh(1) + c*xh(2);
  ap_y = -c * xh_new + d*(1-c)*xh(1) + xh(2);
  xh = [xh_new, xh(1)];
  y(n) = 0.5 * (x(n) - ap_y);  % change to plus for bandreject
end;



function y = aplowpass (x, Wc)
% y = aplowpass (x, Wc)
% Author: M. Holters
% Applies a lowpass filter to the input signal x.
% Wc is the normalized cut-off frequency 0<Wc<1, i.e. 2*fc/fS.
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

c = (tan(pi*Wc/2)-1) / (tan(pi*Wc/2)+1);
xh = 0;
for n=1:length(x)
  xh_new = x(n) - c*xh;
  ap_y = c * xh_new + xh;
  xh = xh_new;
  y(n) = 0.5 * (x(n) + ap_y);  % change to minus for highpass
end;
*/


function y = lowshelving (x, Wc, G)
% y = lowshelving (x, Wc, G)
% Author: M. Holters
% Applies a low-frequency shelving filter to the input signal x.
% Wc is the normalized cut-off frequency 0<Wc<1, i.e. 2*fc/fS
% G is the gain in dB
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

V0 = 10^(G/20); H0 = V0 - 1;
if G >= 0
  c = (tan(pi*Wc/2)-1) / (tan(pi*Wc/2)+1);     % boost
else
  c = (tan(pi*Wc/2)-V0) / (tan(pi*Wc/2)+V0);   % cut
end;
xh = 0;
for n=1:length(x)
  xh_new = x(n) - c*xh;
  ap_y = c * xh_new + xh;
  xh = xh_new;
  y(n) = 0.5 * H0 * (x(n) + ap_y) + x(n);  % change to minus for HS
end;

% Authors: P.Dutilleux, U. Zölzer
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

x=zeros(100,1);x(1)=1; % unit impulse signal of length 100
g=0.5;
b_0=0.5;
b_1=0.5;
a_1=0.7;
xhold=0;yhold=0;
Delayline=zeros(10,1); % memory allocation for length 10
for n=1:length(x);
    yh(n)=b_0*Delayline(10)+b_1*xhold-a_1*yhold; 
    % 1st-order difference equation
    yhold=yh(n);
    xhhold=Delayline(10);
    y(n)=x(n)+g*yh(n);
    Delayline=[y(n);Delayline(1:10-1)];
end;

function y = peakfilt (x, Wc, Wb, G)
% y = peakfilt (x, Wc, Wb, G)
% Author: M. Holters
% Applies a peak filter to the input signal x.
% Wc is the normalized center frequency 0<Wc<1, i.e. 2*fc/fS.
% Wb is the normalized bandwidth 0<Wb<1, i.e. 2*fb/fS.
% G is the gain in dB.
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

V0 = 10^(G/20); H0 = V0 - 1;
if G >= 0
  c = (tan(pi*Wb/2)-1) / (tan(pi*Wb/2)+1);     % boost
else
  c = (tan(pi*Wb/2)-V0) / (tan(pi*Wb/2)+V0);   % cut
end;
d = -cos(pi*Wc);
xh = [0, 0];
for n=1:length(x)
  xh_new = x(n) - d*(1-c)*xh(1) + c*xh(2);
  ap_y = -c * xh_new + d*(1-c)*xh(1) + xh(2);
  xh = [xh_new, xh(1)];
  y(n) = 0.5 * H0 * (x(n) - ap_y) + x(n);
end;
*/


% Authors: P. Dutilleux, U Zölzer
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

x=zeros(100,1);x(1)=1; % unit impulse signal of length 100
BL=0.5;
FB=-0.5;
FF=1;
M=10;
Delayline=zeros(M,1); % memory allocation for length 10
for n=1:length(x);
    xh=x(n)+FB*Delayline(M);
    y(n)=FF*Delayline(M)+BL*xh; 
    Delayline=[xh;Delayline(1:M-1)];
end;

function y=vibrato(x,SAMPLERATE,Modfreq,Width)
% Author: S. Disch
%
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------

ya_alt=0;
Delay=Width; % basic delay of input sample in sec
DELAY=round(Delay*SAMPLERATE); % basic delay in # samples
WIDTH=round(Width*SAMPLERATE); % modulation width in # samples
if WIDTH>DELAY 
  error('delay greater than basic delay !!!');
  return;
end
MODFREQ=Modfreq/SAMPLERATE; % modulation frequency in # samples
LEN=length(x);        % # of samples in WAV-file
L=2+DELAY+WIDTH*2;    % length of the entire delay  
Delayline=zeros(L,1); % memory allocation for delay
y=zeros(size(x));     % memory allocation for output vector
for n=1:(LEN-1)
   M=MODFREQ;
   MOD=sin(M*2*pi*n);
   TAP=1+DELAY+WIDTH*MOD;
   i=floor(TAP);
   frac=TAP-i;
   Delayline=[x(n);Delayline(1:L-1)]; 
   %---Linear Interpolation-----------------------------
   y(n,1)=Delayline(i+1)*frac+Delayline(i)*(1-frac);
   %---Allpass Interpolation------------------------------
   %y(n,1)=(Delayline(i+1)+(1-frac)*Delayline(i)-(1-frac)*ya_alt);  
   %ya_alt=ya(n,1);
   %---Spline Interpolation-------------------------------
   %y(n,1)=Delayline(i+1)*frac^3/6
   %....+Delayline(i)*((1+frac)^3-4*frac^3)/6
   %....+Delayline(i-1)*((2-frac)^3-4*(1-frac)^3)/6
   %....+Delayline(i-2)*(1-frac)^3/6; 
   %3rd-order Spline Interpolation
end  
