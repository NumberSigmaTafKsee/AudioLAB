casino = require('stdsamples_casino')
require('plot')
w = casino.wav_data()
x = casino.load_wave_float("../Data/temp.wav",w)
dft = casino.RDFT32(x:size())
b = casino.VectorXcf(x:size())
casino.dft(dft,x:data(),b:data())
p = plot.Plot_Float()

mag = casino.VectorXf(x:size()/2)
for i=1,x:size()/2 do
    mag[i] = casino.cabsf(b[i])
end
p:setstyle("boxes")
p:plot_x(mag:data(),mag:size(),"audio")
os.execute('sleep 15;')
