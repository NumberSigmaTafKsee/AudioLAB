require('analog_blits_oscillators')
require('stdsamples')
require('plot')
x = analog_blits_oscillators.BlitSaw()
print(x.nHarmonics_)
m = stdsamples.sample_vector(256)
for j=1,10 do
for i=1,256 do
    m[i] = 0.85*x:Tick()
end
end
p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")


x = analog_blits_oscillators.BlitSquare()
m = stdsamples.sample_vector(256)
for j=1,10 do
for i=1,256 do
    m[i] = x:Tick()
end
end
p:plot_x(m:data(),256,"square")

x = analog_blits_oscillators.BlitTriangle()
x:setFrequency(440)
m = stdsamples.sample_vector(256)
for j=1,10 do
for i=1,256 do
    m[i] = x:Tick()    
end
end
p:plot_x(m:data(),256,"triangle")
os.execute("sleep 15;")
