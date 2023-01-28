function x = wfirhp(fc,taps)
  n = floor(taps/2)
  x =  (-n:n)
  x =  sinc(x)-fc*sinc(fc*x)
  x = x.*hamming(taps+1)'
 end
