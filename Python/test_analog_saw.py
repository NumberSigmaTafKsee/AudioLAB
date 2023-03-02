import Analog
import matplotlib.pyplot as plt
import numpy as np

r = Analog.float_vector(256)
r[0] = 100
print(r[0])
x = Analog.BlitSaw()
x.setFrequency(440)
v = np.zeros(256)
for j in range(10):
    for i in range(0,256):
        v[i] = x.Tick()
plt.plot(v)
plt.show()



x = Analog.BlitSquare()
x.setFrequency(440)
v = np.zeros(256)
for j in range(10):
    for i in range(0,256):
        v[i] = x.Tick()/2
plt.plot(v)
plt.show()


svf = Analog.AnalogSVF(44100,500,0.5)
for i in range(0,256):
    v[i] = svf.Tick(v[i])
plt.plot(v)
plt.show()

