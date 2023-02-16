require('stdsamples')
a = stdsamples.complex_vector(16)
b = stdsamples.complex_vector(16)
c = stdsamples.complex_vector(16)
for i=1,16 do
    a[i] = stdsamples.complex(i,0)
end
b:fill(stdsamples.complex(2,0))
fft = stdsamples.FFTPlanComplexFloat(16);
stdsamples.fft(fft,a:data(),b:data())
require('plot')
--stdsamples.ifft(fft,b:data(),c:data());
p = plot.Plot_Float()

p:data_style("boxes")
p:setstyle("boxes")
p:boxwidth(0.8,"relative")
p:set_fill_solid(1.0)

r = stdsamples.sample_vector(16)
for i=1,16 do
    r[i] = stdsamples.cabsf(b[i])
end
p:plot_x(r:data(),r:size(),"fft")
os.execute("sleep 60;")
