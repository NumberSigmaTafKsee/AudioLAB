casino = require('audiodsp_casino')
require('plot')
w = casino.wav_data()
x = casino.load_wave_float("Data/temp.wav",w)
dft = casino.CDFT32(1024)
v = casino.VectorXf(1024)
out = casino.VectorXcf(x:size())
casino.dft(dft,x:data(),out:data())
p = plot.Plot_Float()

mag = casino.VectorXf(x:size()/2)
for i=1,x:size()/2 do
    mag[i] = casino.cabsf(out[i])
end
p:setstyle("boxes")
p:plot_x(mag:data(),mag:size(),"audio")
os.execute('sleep 15;')
