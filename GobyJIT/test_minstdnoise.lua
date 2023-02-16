require('stdnoise')
r = stdnoise.Default()
print("Default Random",r:min(),r:max())
for i=1,10 do io.write(r:generate(),',') end
print()
r = stdnoise.MinStdNoise()
for i=1,10 do io.write(r:generate(),',') end
print()
for i=1,10 do io.write(r:randint(-5,5),',') end
print()
for i=1,10 do io.write(r:rand(),',') end
print()
for i=1,10 do io.write(r:uniform_int_distribution(-10,10),",") end
print()
for i=1,10 do io.write(r:uniform_real_distribution(-2.4,2.4),",") end
print()
for i=1,10 do io.write(r:binomial_distribution(10,0.4),",") end
print()
for i=1,10 do io.write(r:cauchy_distribution(0.78,1.45),",") end
print()
for i=1,10 do io.write(r:chi_squared_distribution(3.0),",") end
print()
