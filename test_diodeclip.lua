-- use xoscope to see the waveform on the device
require('audiosystem')
require('soundwave')
require('vaanalogsvf')
require('vablitoscillator')
require('amplifiers')
require('simpleresampler')
require('functions')
require('vadiodeclipper')

adsr = soundwave.ADSR()
adsr:setAttackRate(0.02*44100)
adsr:setDecayRate(0.3*44100)
adsr:setReleaseRate(0.2*44100)
--sine = soundwave.BandlimitedOscillator(44100,soundwave.BandlimitedOscillator.SQUARE)
saw  = vablitoscillator.BlitSaw()
filt = vaanalogsvf.AnalogSVF(44100,400,0.5)
adsr:gate(0)
freq = 220
rsmp = simpleresampler.SimpleResampler()
rsmp:setup(44100,2)
tmpbuf = simpleresampler.float_vector(512)
A = soundwave.float_vector(256)
X = soundwave.float_vector(256)
Y = soundwave.float_vector(256)
lfo1 = functions.SineGenerator(1.0)
lfo2 = functions.SineGenerator(0.1)
lfo3 = functions.SineGenerator(0.01)
lfo1:setRange(0.2,0.8)
lfo2:setRange(0.1,0.9)
lfo3:setRange(0.3,0.7)
diode = vadiodeclipper.diode_clipper()
diode:setSampleRate(44100)

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
    local buffer = new_buffer(output)    
    saw:setFrequency(freq)
    saw:ProcessBlock(frames,output,output) 
    diode:ProcessBlock(frames,output,output)
    rsmp:up(256,output,tmpbuf:data())  
    for i=1,tmpbuf:size() do tmpbuf[i] = amplifiers.assymetric_sigmoid(tmpbuf[i],5) end
    rsmp:down(256,tmpbuf:data(),output)
    for i=1,256 do
        A[i] = lfo1:Tick()
        X[i] = lfo2:Tick()
        Y[i] = lfo3:Tick()        
    end
    filt:ProcessBlock(frames,output,output,A:data(),X:data(),Y:data())  
    adsr:Process(frames,output,output)   
end 

num_keys=0
function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end
function note_on(c,n,v)      
    local f = math.pow(2.0, (n-69)/12)*440.0            
    freq = f    
    if(num_keys == 0) then
        adsr:gate(1)    
    end
    num_keys = num_keys+1
end
function note_off(c,n,v) 
    num_keys = num_keys-1
    if(num_keys <= 0) then
        adsr:gate(0)    
        num_keys = 0;
    end       
end
function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 102) then
        filt:setCutoff(22050.0*d2/127)
    elseif(d1 == 103) then
        q = 10*d2/127
        if(q >= 0.5) then
            filt:setQ(q)
        end
    end
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
