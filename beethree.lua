-- use xoscope to see the waveform on the device
require('audiosystem')
require('soundwave')
require('stk')
require('lv2plugin')
require('faustfx')

stk.setRawwavePath("Data/rawwaves")
stk.setSampleRate(44100.0)
b3 = stk.BeeThree()
freq = 440.0
lv2plugins = lv2plugin.LV2Plugins()
gverb = lv2plugins:LoadPlugin("http://plugin.org.uk/swh-plugins/gverb")
flanger = faustfx.FaustFX("flanger.dsp")

function new_buffer(p)
    local b = {}
    b.buffer = p
    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

function noise(input,output,frames)            
    local outbuf = new_buffer(output)    
    b3:setFrequency(freq)    
    for i=0,frames-1 do
        local out = b3:tick(0)        
        outbuf[i] = out                
    end
    flanger:Run(frames,output,output)
    gverb:Run(frames,output,output)
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
    b3:noteOn(f,v/127)    
    numpress = numpress+1
end
function note_off(c,n,v)        
    numpress = numpress-1
    if(numpress <= 0) then 
        numpress = 0;
        b3:keyOff()    
    end
end

-- app callback, midi handling and logic
-- isAudioRunning shuld be atomic
-- either block audio or wait until finished
-- this is run every 10ms, or can be changed in portaudio.i
function callback()
    print('hi')
end 

function randomize()        
    b3:setRatio(0,math.random()*4);
    b3:setRatio(1,math.random()*4);
    b3:setRatio(2,math.random()*4);
    b3:setRatio(3,math.random()*4);
    b3:setGain(0,math.random()*4);
    b3:setGain(1,math.random()*4);
    b3:setGain(2,math.random()*4);
    b3:setGain(3,math.random()*4);
    b3:setModulationSpeed(math.random()*10);
    b3:setModulationDepth(math.random());
    --b3:setControl1(math.rand()*10);
    --b3:setControl2(math.rand()*10);
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
audiosystem.InitAudioDevice(device,-1,1,44100,256)
audiosystem.RunAudio()
audiosystem.StopAudio()
