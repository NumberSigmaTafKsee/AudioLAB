require('kfr2')

v = kfr2.SampleVector(10)
for i=1,10 do v[i] = i end
v:print()
x = kfr2.acorr(v)
x:print()
x = kfr2.xcorr(v,x)
x:print()