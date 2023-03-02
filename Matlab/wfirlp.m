function x = wfirlp(fc,taps)
  n = floor(taps/2)
  x =  fc * sinc(fc*(-n:n))
  x = x.*hamming(taps+1)'
 end
