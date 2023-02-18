require('spucecheby2bp')
require('sndfile')
s = sndfile.SndFileReaderFloat("baby_elephant.wav")
print(s:channels())
v = sndfile.float_vector(s:frames()*s:channels())
s:read(v:size(),v)
bp = spucecheby2bp.BandpassChebyshev2(4,350,450,10,s:samplerate())
for i=1,v:size() do v[i] = bp:Tick(v[i]) end
o = sndfile.SndFileWriterFloat("test.wav",0x10006,s:channels(),s:samplerate())
o:write(v)