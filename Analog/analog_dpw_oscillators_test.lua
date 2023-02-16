require('analog_dpw_oscillators')
require('stdsamples')
require('plot')

x = analog_dpw_oscillators.DPWTriangle()
x:setFrequency(440)
m = stdsamples.sample_vector(256)

for j=1,10 do
for i=1,256 do
    m[i] = x:Tick()      
end
end

print('done')

p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")
os.execute('sleep 15;')