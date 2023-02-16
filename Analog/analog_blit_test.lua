require('Analog')
require('stdsamples')
require('plot')

x = Analog.blitSaw()
m = stdsamples.sample_vector(256)
for j=1,10 do
for i=1,256 do
    m[i] = x:Tick()
end
end
p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")


x = Analog.blitSquare()
m = stdsamples.sample_vector(256)
for j=1,10 do
for i=1,256 do
    m[i] = x:Tick()
end
end
p:plot_x(m:data(),256,"square")

x = Analog.blitTriangle()
x:setFrequency(440)
m = stdsamples.sample_vector(256)
for j=1,10 do
for i=1,256 do
    m[i] = x:Tick()    
end
end
p:plot_x(m:data(),256,"triangle")
os.execute("sleep 15;")
