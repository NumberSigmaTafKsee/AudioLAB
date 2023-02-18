require('audiodsp_adsr')
require('plot')
require('stdsamples')
require('stdsamples_function_generator')

p = plot.Plot_Float()
p:setstyle("lines")
v = stdsamples.sample_vector_float(256)
a = audiodsp_adsr.ADSR(0.01,0.1,0.9,0.01)
a:noteOn()
for i=1,128 do
	v[i] = a:Tick()
end
a:noteOff()
for i=128,256 do
	v[i] = a:Tick()
end
p:plot_x(v:data(),v:size(),"amp")
os.execute("sleep 15;")

