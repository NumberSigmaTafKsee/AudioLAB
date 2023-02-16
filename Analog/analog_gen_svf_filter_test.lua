require('analog_dpw_oscillators')
require('analog_gen_svf_filter')
require('stdsamples')
require('plot')

x = analog_dpw_oscillators.DPWTriangle()
x:setFrequency(440)
f = analog_gen_svf_filter.GeneralSVF()
f:calcCoefs(0.5,0,5000.0/44100.0)
m = stdsamples.sample_vector(256)

for j=1,10 do
for i=1,256 do
    m[i] = x:Tick()    
    m[i] = f:process(m[i])  
end
end

print('done')

p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")
os.execute('sleep 15;')