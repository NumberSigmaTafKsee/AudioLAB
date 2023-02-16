require('kfr2')
fmt = kfr2.audio_format()
wav = kfr2.wav_load('Data/fairytale.wav',fmt)
bqs = kfr2.biquad_params_vector(1)
ap  = kfr2.biquad_lowpass(500.0/44100.0,2)
bqs[1] = ap
output = kfr2.biquad_filter(bqs)
output:apply(wav)
kfr2.wav_save(wav,"test.wav",fmt.channels, fmt.samplerate)
