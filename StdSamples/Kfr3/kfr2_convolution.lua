require('kfr2')

fmt = kfr2.audio_format()
wav = kfr2.wav_load('Data/redwheel.wav',fmt)
print(fmt.channels, fmt.samplerate)
h   = kfr2.wav_load("01 Halls 1 Large Hall.wav",fmt)
m   = kfr2.deinterleave(h)
L   = m[1]
R   = m[2]
c   = kfr2.ConvolutionFilter(L)
out = c:Process(wav)
kfr2.wav_save(out,"reverb.wav",1,44100)