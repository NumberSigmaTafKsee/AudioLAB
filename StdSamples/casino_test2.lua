casino = require('stdsamples_casino')
require('plot')
w = casino.wav_data()
x = casino.load_wave_float("../Data/temp.wav",w)
n = 2^math.ceil(math.log(x:size())/math.log(2))
print(n)
dft = casino.RFFT32(n)
v = casino.VectorXf(n)
out = casino.VectorXcf(n)
casino.fft(dft,x:data(),out:data())
p = plot.Plot_Float()

mag = casino.VectorXf(n/2)
for i=1,n/2 do
    mag[i] = casino.cabsf(out[i])
end
p:setstyle("boxes")
p:plot_x(mag:data(),mag:size(),"audio")
os.execute('sleep 15;')
