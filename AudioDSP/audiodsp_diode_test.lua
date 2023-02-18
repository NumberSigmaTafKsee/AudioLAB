require('audiodsp_blits_oscillators')
require('audiodsp_chamberlin_svf_filter')
require('audiodsp_diode')
require('stdsamples')
require('plot')
x = audiodsp_blits_oscillators.BlitSaw()
svf = audiodsp_chamberlin_svf_filter.ChamberlinSVF(44100,5000,0.5)
svf:setCutoff(1000);
m = stdsamples.sample_vector(256)

for j=1,10 do
for i=1,256 do
    m[i] = 0.9*x:Tick()    
    m[i] = audiodsp_diode.Diode(m[i])
    local x = svf:Tick(m[i])
    m[i] = math.sqrt(math.sqrt(temp*temp-x*x))
    temp = x
end
end
max = stdsamples.maxf(

p = plot.Plot_Float()
p:setstyle("lines")
p:plot_x(m:data(),256,"saw")
os.execute('sleep 15;')
