require('vf')

a = vf.Vector(100)
c = vf.Vector(100)
d = vf.Vector(100)
c:fill(1.0)
d:fill(1/1.25)
for i=0,99 do a[i] = i/100  end
a:upload_device()
r = c + vf.exp(a*a)/2
r:print()

