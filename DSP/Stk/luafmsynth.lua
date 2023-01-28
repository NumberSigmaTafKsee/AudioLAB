-- use xoscope to see the waveform on the device
require('audiosystem')
require('soundwave')
require('fmsynth')
require('lv2plugin')
require('faustfx')


fm = fmsynth.FMSynth(44100,64)
preset = fmsynth.FMPreset("presets/old_school_organ.fmp")
fm:load_preset(preset)
freq = 440.0
lv2plugins = lv2plugin.LV2Plugins()
gverbL = lv2plugins:LoadPlugin("http://plugin.org.uk/swh-plugins/gverb")
flangerL = faustfx.FaustFX("flanger.dsp")
gverbR = lv2plugins:LoadPlugin("http://plugin.org.uk/swh-plugins/gverb")
flangerR = faustfx.FaustFX("flanger.dsp")

print('go')
function new_buffer(p)
    local b = {}
    b.buffer = p
    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

sampleRate=44100
nframes=256
L = fmsynth.float_vector(nframes)   
R = fmsynth.float_vector(nframes)   

function noise(input,output,frames)            
    local outbuf = new_buffer(output) 
    for i=1,nframes do 
        L[i] = 0
        R[i] = 0
    end
    fm:render(L:data(),R:data(),frames)    
    flangerL:Run(frames,L:data(),L:data())
    gverbL:Run(frames,L:data(),L:data())
    flangerR:Run(frames,R:data(),R:data())
    gverbR:Run(frames,R:data(),R:data())
    for i=0,frames-1 do        
        outbuf[i*2] = L[i+1]
        outbuf[i*2+1]=R[i+1]
    end
end 

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end
numpress = 0
function note_on(c,n,v)    
    local f = math.pow(2.0, (n-69)/12)*440.0            
    freq = f
    fm:note_on(n,v)    
    numpress = numpress+1
end
function note_off(c,n,v)        
    fm:note_off(n)    
end

-- app callback, midi handling and logic
-- isAudioRunning shuld be atomic
-- either block audio or wait until finished
-- this is run every 10ms, or can be changed in portaudio.i
function callback()
    print('hi')
end 


audiosystem.Init()
audiosystem.Pm_Initialize()

audiosystem.set_note_on_func(note_on)
audiosystem.set_note_off_func(note_off)

for i=0,audiosystem.GetNumMidiDevices()-1 do 
    print(i,audiosystem.GetMidiDeviceName(i))
end

audiosystem.set_audio_func(noise)
--audiosystem.set_callback_func(callback)
device=14
audiosystem.Pa_Initialize()
for i=0,audiosystem.GetNumAudioDevices()-1 do 
    if( audiosystem.GetAudioDeviceName(i) == 'jack') then        
        device = i 
        goto done
    end    
end
::done::
audiosystem.InitMidiDevice(1,3,3)
audiosystem.InitAudioDevice(device,-1,2,sampleRate,nframes)
audiosystem.RunAudio()
audiosystem.StopAudio()
