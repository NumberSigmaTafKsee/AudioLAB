function x = wfirls(G,fc,taps)
  n = floor(taps/2)
  lp = wfirlp(fc,taps)
  hp = wfirhp(fc,taps)
  x = G*lp + hp
  x = x.*hamming(taps+1)'
end
