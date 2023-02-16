require('audiofft')
require('sndfile')

read = sndfile.SoundFileReader()
read:open("test.wav")
print(read:frames())

fft = audiofft.AudioFFT()
fft:init(1024)
v = sndfile.float_vector(1024)
real = sndfile.float_vector(1024)
imag = sndfile.float_vector(1024)
read:read(1024,v)
fft:fft(v:data(),real:data(),imag:data())

mag = 0
for i=0,511 do
    local x = math.sqrt(real[i]*real[i]+imag[i]*imag[i])
    if(x > mag) then mag = x end 
end 
print(mag/512)