require('DspFilters')
require('sndfile')

s = sndfile.SndFileReaderFloat("baby_elephant.wav")
v = sndfile.float_vector(s:size())
s:read(v:size(),v:data())
f = DspFilters.RBJLowPass(1024)
p = DspFilters.Params();
p[0] = s:samplerate();
p[1] = 1000;
p[2] = 1.25;
f:setParams(p)
f:ProcessBlock(v:size(),v:data(),v:data())
o = sndfile.SndFileWriterFloat("test.wav",0x10006,s:channels(),s:samplerate())
o:write(v)