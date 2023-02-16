require('analog_blits_oscillators')
require('analog_chamberlin_svf_filter')
require('analog_diode_clipper')
require('stdsamples')
require('plot')
x = analog_blits_oscillators.BlitSaw()
svf = analog_chamberlin_svf_filter.ChamberlinSVF(44100,1000,0.5)
d = analog_diode_clipper.DiodeClipper()
m = stdsamples.sample_vector(256)
for j=1,10 do
for i=1,256 do
    m[i] = 1.5*x:Tick()
    m[i] = d:Tick(m[i])
    m[i] = svf:Tick(m[i])
end
end

p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")
os.execute('sleep 15;')