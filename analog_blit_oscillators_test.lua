require('analog_blit_oscillators')
require('stdsamples')
require('plot')
x = analog_blit_oscillators.blitSaw()
m = stdsamples.sample_vector(256)
for j=1,10 do
    x:ProcessSIMD(256,m:data())
end
p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")


x = analog_blit_oscillators.blitSquare()
m = stdsamples.sample_vector(256)
for j=1,10 do
    x:ProcessSIMD(256,m:data())
end
p:plot_x(m:data(),256,"square")

x = analog_blit_oscillators.blitTriangle()
x:setFrequency(4000)
m = stdsamples.sample_vector(256)
for j=1,10 do
    x:ProcessSIMD(256,m:data())
end
p:plot_x(m:data(),256,"triangle")
os.execute("sleep 15;")
