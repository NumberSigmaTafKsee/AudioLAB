function x = wfirbs(fl,fu,taps)
  n = floor(taps/2)
  lp= wfirlp(fl,taps)
  hp= wfirhp(fu,taps)
  x = lp+hp
  x = x.*hamming(taps+1)'
end
