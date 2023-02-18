require('audiodsp_amplifier')
require('plot')
require('stdsamples')
require('stdsamples_function_generator')

p = plot.Plot_Float()
p:setstyle("lines")
v = stdsamples.generate_sinf(440,44100,256)
a = audiodsp_amplifier.Amplifier()
a:ProcessSIMD(v:size(),v:data(),v:data())
p:plot_x(v:data(),v:size(),"amp")
os.execute("sleep 15;")

