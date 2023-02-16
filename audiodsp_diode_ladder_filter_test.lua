require('analog_blits_oscillators')
require('analog_diode_ladder_filter')
require('analog_diode_clipper')
require('stdsamples')
require('plot')

x = analog_blits_oscillators.BlitSaw()
f = analog_diode_ladder_filter.DiodeLadderFilter(44100)
f:setCutoff(5)
f:setResonance(0.5)
d = analog_diode_clipper.DiodeClipper()
m = stdsamples.sample_vector(256)

for i=1,256 do
    m[i] = x:Tick()    
    m[i] = 2*f:Tick(m[i])    
end
print('done')

p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")
os.execute('sleep 15;')