% EQ types: lp, hp, bp, ap, highShelf, lowShelf, peak, notch
% lp -- lowpass
% hp -- highpass
% bp -- bandpass
% ap -- allpass
function [arr] = RbjcookEQ(Fs,f0,dBgain,Q,filterType)


if (strcmp(filterType,'peak') || strcmp(filterType,'lowShelf') || strcmp(filterType,'highShelf'))
  A  = sqrt(10^(dBgain/40));
else
  A  = sqrt(10^(dBgain/20));
end

w0 = 2*pi*f0/Fs;
alpha = sin(w0)/(2*Q);


switch filterType
  case 'lp'
    b0 =  (1 - cos(w0))/2;
    b1 =   1 - cos(w0);
    b2 =  (1 - cos(w0))/2;
    a0 =   1 + alpha;
    a1 =  -2*cos(w0);
    a2 =   1 - alpha;
  case 'hp'
    b0 =  (1 + cos(w0))/2;
    b1 = -(1 + cos(w0));
    b2 =  (1 + cos(w0))/2;
    a0 =   1 + alpha;
    a1 =  -2*cos(w0);
    a2 =   1 - alpha;
  case 'bp'
    b0 =   alpha;
    b1 =   0;
    b2 =  -alpha;
    a0 =   1 + alpha;
    a1 =  -2*cos(w0);
    a2 =   1 - alpha;
  case 'notch'
    b0 =   1;
    b1 =  -2*cos(w0);
    b2 =   1;
    a0 =   1 + alpha;
    a1 =  -2*cos(w0);
    a2 =   1 - alpha;
  case 'ap'
    b0 =   1 - alpha;
    b1 =  -2*cos(w0);
    b2 =   1 + alpha;
    a0 =   1 + alpha;
    a1 =  -2*cos(w0);
    a2 =   1 - alpha;
  case 'peak'
    b0 =   1 + alpha*A;
    b1 =  -2*cos(w0);
    b2 =   1 - alpha*A;
    a0 =   1 + alpha/A;
    a1 =  -2*cos(w0);
    a2 =   1 - alpha/A;
  case 'lowShelf'
    b0 =    A*( (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha );
    b1 =  2*A*( (A-1) - (A+1)*cos(w0)                   );
    b2 =    A*( (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha );
    a0 =        (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha;
    a1 =   -2*( (A-1) + (A+1)*cos(w0)                   );
    a2 =        (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha;
  case 'highShelf'
    b0 =    A*( (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha );
    b1 = -2*A*( (A-1) + (A+1)*cos(w0)                   );
    b2 =    A*( (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha );
    a0 =        (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha;
    a1 =    2*( (A-1) - (A+1)*cos(w0)                   );
    a2 =        (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha;
end

c0 = b0/a0;
c1 = b1/a0;
c2 = b2/a0;
c3 = a1/a0;
c4 = a2/a0;

arr = [c0 c1 c2 1 c3 c4];

