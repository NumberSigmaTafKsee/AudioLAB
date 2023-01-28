function x = wfirhs(G,fc,taps)
  n = floor(taps/2)
  lp= wfirlp(fc,taps)
  hp= wfirhp(fc,taps)
  x = lp+G*hp
  x = x.*hamming(taps+1)'
end
