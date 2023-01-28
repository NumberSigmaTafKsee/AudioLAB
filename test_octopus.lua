require('octopus')
interp = octopus.OctaveInterpreter()
v = octopus.OctopusColVectorXf(3)

v:fill(10)
octopus.display(v)

c = v + 10
octopus.display(c)

c = v - 10
octopus.display(c)

c = v * 10
octopus.display(c)

c = v / 10
octopus.display(c)

c = v + v
octopus.display(c)

c = v - v
octopus.display(c)

c = v * v
octopus.display(c)

c = v / v
octopus.display(c)