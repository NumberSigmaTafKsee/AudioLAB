require('octopus')
require('stdsamples')
i = octopus.OctaveInterpreter()
v = octopus.OctopusRowVectorXf(10)
for i=1,10 do v[i] = i end
s = stdsamples.sample_vector(10)
s:copy(v:size(),v:data())
s:print()
l = octopus.OctopusValueList()
l[0] = 0
l[1] = 1
x = octopus.octave_linspace(l,1)
r = x[0]:getFloatRowVector() 
r:print()
r = r * 2 * math.pi 

--l = octopus.OctopusValueList()
--l[0] = r
m = octopus.octave_sin(r)
--m = l[0]:getFloatRowVector()
m:print()