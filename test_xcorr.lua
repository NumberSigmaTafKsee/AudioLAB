require('kfr2')

v = kfr2.SampleVector(10)
for i=1,10 do v[i] = i end
x = kfr2.xcorr(v,v)
x:print()