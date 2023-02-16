require('soundwave')
require('audiosystem')

osc = soundwave.BandlimitedOscillator(44100,soundwave.SINE)

function new_buffer(p)
    local b = {}
    b.buffer = p

    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

function sound(input,output,frames)
    local out = new_buffer(output)
    for i=0,frames-1 do 
        local sample = osc:tick()
        out[i*2] = sample
        out[i*2+1] = sample        
    end 
end

audiosystem.set_audio_func(sound)
--luapa.set_callback_func(callback)
audiosystem.Pa_Initialize()
for i=0,audiosystem.GetNumAudioDevices()-1 do 
    print(i,audiosystem.GetAudioDeviceName(i))
end
-- use jack on my system it is device 6
-- use pulse on my system it is device 14
-- no input is used (-1)
audiosystem.InitAudioDevice(10,-1,2,44100,256)
audiosystem.RunAudio()
audiosystem.Sleep(1000)
audiosystem.StopAudio()