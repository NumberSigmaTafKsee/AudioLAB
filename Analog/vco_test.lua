require('VCO')
require('stdsamples')
require('plot')
x = VCO.VCO(44100.0,VCO.VCO.TRAPEZOID_VARIABLE)
m = stdsamples.sample_vector(256)
for j=1,10 do
    x:ProcessSIMD(256,m:data())
end
p = plot.Plot_Double()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")
os.execute("sleep 15;")
