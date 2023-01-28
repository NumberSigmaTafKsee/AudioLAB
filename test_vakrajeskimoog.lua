-- use xoscope to see the waveform on the device
require('audiosystem')
require('vector')
require('adsr')
require('vaanalogsvf')
require('vablitoscillator')
require('amplifiers')
require('simpleresampler')
require('functions')
require('vadiodeclipper')
require('vakrajeskimoog')

adsr1 = adsr.ADSR()
adsr1:setAttackRate(0.02*44100)
adsr1:setDecayRate(0.3*44100)
adsr1:setReleaseRate(0.2*44100)


adsr2 = adsr.ADSR()
adsr2:setAttackRate(1.5*44100)
adsr2:setDecayRate(0.5*44100)
adsr2:setReleaseRate(0.5*44100)
adsr2:setSustainLevel(0.5)
adsr2.max = 2

saw  = vablitoscillator.BlitSaw()
filt = vakrajeskimoog.KrajeskiMoog(44100)
rsmp = simpleresampler.SimpleResampler()
rsmp:setup(44100,2)
tmpbuf = simpleresampler.float_vector(512)
A = vector.float_vector(256)
X = vector.float_vector(256)
Y = vector.float_vector(256)
lfo1 = functions.SineGenerator(1.0)
lfo2 = functions.SineGenerator(0.1)
lfo3 = functions.SineGenerator(0.01)
lfo1:setRange(0.2,0.8)
lfo2:setRange(0.1,0.9)
lfo3:setRange(0.3,0.7)
diode = vadiodeclipper.diode_clipper()
diode:setSampleRate(44100)
diode:initialise()

function new_buffer(p)
    local b = {}
    b.buffer = p
    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

freq=400
fc = 440
qr = 0.5
function noise(input,output,frames)        
    local buffer = new_buffer(output)    
    saw:setFrequency(freq)    
    filt:SetCutoff(fc)
    filt:SetResonance(4*qr)
    saw:ProcessBlock(frames,output,output) 
    diode:ProcessBlock(frames,output,output)
    --rsmp:up(256,output,tmpbuf:data())  
    --for i=1,tmpbuf:size() do tmpbuf[i] = amplifiers.assymetric_sigmoid(tmpbuf[i],5) end
    --rsmp:down(256,tmpbuf:data(),output)    
    for i=1,frames do 
        A[i] = adsr2:Tick()
        X[i] = lfo1:Tick()
        Y[i] = lfo2:Tick()
        buffer[i-1] = adsr1:Tick()*filt:Tick(buffer[i-1],A[i],X[i],Y[i])
    end
    --filt:ProcessBlock(frames,output,output,A:data())
    --adsr1:ProcessBlock(frames,output,output)   
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
        adsr1:noteOn() 
        adsr2:noteOn()    
    end
    num_keys = num_keys+1
end
function note_off(c,n,v) 
    num_keys = num_keys-1
    if(num_keys <= 0) then
        adsr1:noteOff()
        adsr2:noteOff()   
        num_keys = 0;
    end       
end
function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 102) then
        fc = (22050/2.0*d2/127)
    elseif(d1 == 103) then
        qr = d2/127                
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
