require('analog_dpw_oscillators')
require('analog_improved_moog_filter')
require('stdsamples')
require('plot')

x = analog_dpw_oscillators.DPWSaw()
x:setFrequency(440)
f = analog_improved_moog_filter.ImprovedMoog(44100.0)
m = stdsamples.sample_vector(256)

for j=1,10 do
for i=1,256 do
    m[i] = x:Tick()    
    m[i] = f:Tick(m[i])  
end
end

print('done')

p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")
os.execute('sleep 15;')