require('kfr2')
require('kfr_filters')
fmt = kfr2.audio_format()
wav = kfr2.wav_load('Data/fairytale.wav',fmt)
print(fmt.samplerate)
lp  = kfr_filters.BiquadLowpassFilter(4,400.0,fmt.samplerate,0.707)
out = kfr2.univector(wav:size())
lp:ProcessBlock(wav:size(),wav:data(),out:data())
kfr2.wav_save(out,"test.wav",fmt.channels, fmt.samplerate)