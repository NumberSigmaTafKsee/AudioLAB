require('sndfile')
require('gist')

s = sndfile.SndFileReaderFloat("RhythmGuitar.wav")
o = sndfile.SndFileWriterFloat('testsynth.wav',0x10006,s:channels(),s:samplerate())
v = sndfile.float_vector(s:frames()*s:channels())
s:read(v:size(),v)
-- non-rectangular windows will distort it
g = gist.GistFloat(1024,s:samplerate(),gist.RectangularWindow)
x = sndfile.float_vector(1024)
for i=1,v:size()/1024 do    
    for j=1,1024 do x[j] = v[i*1024+j] end    
    g:processAudioFrame(x)
    m = g:getMagnitudeSpectrum()    
    y = g:synthesize()       
    o:write(y)
end