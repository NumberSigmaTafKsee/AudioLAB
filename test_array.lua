require('se')

a1 = se.FloatArray(100)
a1:random() 
a1:print()
a2 = se.FloatArray(100)
a2:fill(10)
r = a1 * a2 
r:print()
