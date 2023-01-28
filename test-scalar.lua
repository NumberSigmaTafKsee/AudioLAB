require('se')

function const(a) return se.FloatScalar(a) end
s1 = se.FloatScalar(10)
s2 = se.FloatScalar(20)
r = s1 + s2 
print(r.val)
s1:set_value(2*math.pi)
r = s1:cos() 
print(r.val)
r = s1:sin() 
print(r.val)
r = s1:tan() 
print(r.val)
s1:set_value(math.pi/2)
r = s1:cos() 
print(r.val)
r = s1:sin() 
print(r.val)
r = s1:tan() 
print(r.val)
r = const(-10)*s1*s2 + 1000
print(r.val)
