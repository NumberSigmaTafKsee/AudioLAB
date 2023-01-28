function [hmp] = minphaseir(h) 
%
% MINPHASEIR - Convert a real impulse response to its 
%              minimum phase counterpart
%
% USAGE:
%
%       [hmp] = minphaseir(h) 
%
% where
%
% h   = impulse response (any length - will be zero-padded)
% hmp = min-phase impulse response (at zero-padded length)

nh = length(h);
nfft = 2^nextpow2(5*nh);
Hzp = fft(h,nfft);
Hmpzp = exp( fft( fold( ifft( log( clipdb(Hzp,-100) )))));
hmpzp = ifft(Hmpzp);
hmp = hmpzp(1:nh);
relerr = norm(imag(hmp))/norm(hmp);
if relerr > 0.01
 disp(sprintf(...
  'minphaseir: WARNING: norm(imag(hmp))/norm(hmp) = %f',relerr));
end
hmp = real(hmp);
