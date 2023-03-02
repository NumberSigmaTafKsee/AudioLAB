clc
kp = 3;          % passband ripple
ks = 60;          % stop attenuation
fp = 40;          %passband frequency
fs = 1000;        %sampling frequency
Fs = 150;         % stop band frequency
wp = fp/(fs/2);
ws = Fs/(fs/2);
[N wc] = buttord(wp, ws, kp, ks);
disp([N wc])
