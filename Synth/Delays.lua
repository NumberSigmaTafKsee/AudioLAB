-- use xoscope to see the waveform on the device
require('audiosystem')
--require('soundwave')
require('stk')
require('lv2plugin')
require('faustfx')
require('canalog')
require('fxdelays')

stk.setRawwavePath("rawwaves")
stk.setSampleRate(44100.0)
osc  = canalog.CQBLimitedOscillator()
osc.m_dOscFo = 440.0
osc.m_uWaveform = canalog.COscillator.SAW1
osc:setSampleRate(44100)
osc:startOscillator()
filt = canalog.CKThreeFiveFilter()
filt.m_uFilterType = canalog.CFilter.LPF2
filt:setFcMod(1000)
filt:setQControl(0.5)
filt:setSampleRate(44100)
adsr = stk.ADSR()
adsr:setAttackRate(0)
adsr:setAllTimes(0.1,0.2,0.7,0.2)
lv2plugins = lv2plugin.LV2Plugins()
gverb = lv2plugins:LoadPlugin("http://plugin.org.uk/swh-plugins/gverb")
flanger = faustfx.FaustFX("flanger.dsp")
delay   = fxdelays.AudioDelay()
params  = fxdelays.AudioDelayParameters()
params.leftDelay_mSec  = 400.0
params.rightDelay_mSec = 410.0
params.feedback_Pct = 80.0
delay:createDelayBuffers(44100,2000)
delay:setParameters(params)

freq = 440.0
numpress = 0
fc = 440.0
q  = 0.5

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
    osc.m_dOscFo   = freq
    osc:update()
    
    local f = math.pow(127.0,(fc-1.0));
    filt.m_dFcControl = f*22050
    filt.m_dQControl  = 9.9*q    
    filt:setQControl(10*q)
    filt:setFcMod(0)
    filt:update()
    for i=0,frames-1 do        
        local out = osc:doOscillate()   
        out = filt:doFilter(out)     
        out = out * adsr:tick()                 
        outbuf[i] = out                
    end
    delay:processMono(frames,output,output)    
    flanger:Run(frames,output,output)    
    gverb:ProcessBlock(frames,output,output)
end 

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end

function note_on(c,n,v)    
    local f = math.pow(2.0, (n-69)/12)*440.0            
    freq = f    
    adsr:keyOn()    
    numpress = numpress+1
end
function note_off(c,n,v)        
    numpress = numpress-1
    if(numpress <= 0) then 
        numpress = 0;
        adsr:keyOff()
    end
end

function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 102) then        
        fc = d2/127                
    elseif(d1 == 103) then
        q = d2/127                
    end
end


audiosystem.Init()
audiosystem.Pm_Initialize()

audiosystem.set_note_on_func(note_on)
audiosystem.set_note_off_func(note_off)
audiosystem.set_control_change_func(control)


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
