require('audiodsp_blit_oscillators')
require('stdsamples')
require('plot')
x = audiodsp_blit_oscillators.BlitSaw()
m = stdsamples.sample_vector(256)
for j=1,10 do
    x:ProcessSIMD(256,nil,m:data())
end
p = plot.Plot_Float()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")


x = audiodsp_blit_oscillators.BlitSquare()
m = stdsamples.sample_vector(256)
for j=1,10 do
    x:ProcessSIMD(256,nil,m:data())
end
p:plot_x(m:data(),256,"square")

x = audiodsp_blit_oscillators.BlitTriangle()
x:setFrequency(4000)
m = stdsamples.sample_vector(256)
for j=1,10 do
    x:ProcessSIMD(256,nil,m:data())
end
p:plot_x(m:data(),256,"triangle")
os.execute("sleep 15;")
