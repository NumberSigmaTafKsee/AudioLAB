import PyAnalogBlitOscillators
import PyPlot
import time

x = PyAnalogBlitOscillators.Blit2SawOscillator_f32()
v = PyAnalogBlitOscillators.float_vector(256)
for i in range(1,256):
    v[i] = x.Tick()
p = PyPlot.Plot_Float()
p.setstyle("lines")
p.plot_x(v.data(),v.size(),"saw")
time.sleep(16)