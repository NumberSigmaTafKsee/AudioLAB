require('luapa')
require('waveTableOsc')
require('sineosc')
require('waves')

v = waveTableOsc.vectorf(4096)
phase = 0
f = 440.0/44100.0


w = waveTableOsc.WaveTableOsc()
pulse_wave(440.0,v,0.99)
w:setFrequency(f)
w:setPhaseOffset(0.0)
w:addWaveTable(4096,v,440.0)
s = SineOsc(44100,4096)

function noise(input,output,frames)
    local a = 0
    for i = 0,2*frames-1,2 do
        w:updatePhase()
        local x = w:getOutput()
        x = s:Generate(220.0)        
        luapa.float_set(output,i,x)
        luapa.float_set(output,i+1,x)        
    end    
end 


luapa.set_audio_func(noise)
luapa.InitAudio(44100,64)
luapa.RunAudio()
luapa.StopAudio()