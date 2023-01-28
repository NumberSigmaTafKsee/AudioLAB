function x = wfirbp(fl,fu,taps)
  n = floor(taps/2)
  x =  (-n:n)
  x2 = fu*sinc(fu*x)
  x1 = fl*sinc(fl*x)
  x = x2 - x1
  x = x.*hamming(taps+1)'
 end
