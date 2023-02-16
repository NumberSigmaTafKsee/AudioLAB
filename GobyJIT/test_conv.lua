require('kfr2')

v = kfr2.SampleVector(32)
for i=1,32 do v[i] = i end
x = kfr2.conv(v,v)
x:print()