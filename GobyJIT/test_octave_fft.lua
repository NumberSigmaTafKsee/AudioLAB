require('octopus')
interp = octopus.OctaveInterpreter()
x = octopus.OctopusRowVectorXd(128)
for i=1,128 do
    x[i] = i
end
l = octopus.OctopusValueList()
l[0] = x
c = octopus.octave_fft(l,1)
m = c[0]:getComplexMatrix()
m:print()
l[0] = m
l = octopus.octave_ifft(l,1)
m = l[0]:getRowVector()
m:print()