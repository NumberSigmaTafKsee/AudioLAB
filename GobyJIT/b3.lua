require('stk')
require('AudioTK')
require('Amplifiers')
require('audio_lfo')

lfo = audio_lfo.LFO(44100)
lfo:setRate(0.00005)
stk.Stk.setRawwavePath("Data/rawwaves")
stk.Stk.setSampleRate(44100)
adsr = stk.ADSR()
chorus = stk.Chorus()
chorus:setModFrequency(12)
diode = AudioTK.MonoHalfTanhShaper()
delay = AudioTK.MonoFDNDelay(44100)
-- it needs to be oversampled to reduce the aliasing
clip  = Amplifiers.MorphClipper(Amplifiers.ClipFunction.SERPENT_CURVE,Amplifiers.ClipFunction.ERFMOIDER)
diode:setPort(AudioTK.MonoHalfTanhShaper.PORT_COEFF,10)
delay:setPort(AudioTK.MonoFDNDelay.PORT_INGAIN,1)
delay:setPort(AudioTK.MonoFDNDelay.PORT_OUTGAIN,1)
delay:setPort(AudioTK.MonoFDNDelay.PORT_FEEDBACK,0.9)
delay:setPort(AudioTK.MonoFDNDelay.PORT_DELAY,0.5*44100)
adsr:setAllTimes(0.1,0.3,0.8,0.2)
b3   = stk.BeeThree()
envs   = {}
voices = {}
inuse  = {}
freqs  = {}
vels   = {}
voicer = stk.Voicer()
for i=1,8 do
    envs[i] = stk.ADSR()
    envs[i]:setAllTimes(0.1,0.3,0.8,0.2)
    voices[i] = stk.BeeThree()
    voices[i]:noteOff(0)
    inuse[i] = false
    freqs[i] = 0
    vels[i] = 0
    voicer:addInstrument(voices[i])
end


fc = 440.0
q  = 0.5

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end

function new_buffer(p)
    local b = {}
    b.buffer = p
    local mt = {} 
    mt.__index = function(b,i) return audiosystem.float_get(b.buffer,i) end
    mt.__newindex = function(b,i,v) audiosystem.float_set(b.buffer,i,v) end 
    setmetatable(b,mt)
    return b
end 

v = AudioTK.double_vector(256)
v2 = AudioTK.double_vector(256)
function noise(input,output,frames)            
    local outbuf = new_buffer(output)         
    --[[
    local total = 0
    for i=1,8 do 
        if(inuse[i] == true)  then total = total + 1 end      
    end
    for i=0,frames-1 do        
        outbuf[i] = 0
        for j=1,8 do
            outbuf[i] = outbuf[i] + envs[j]:tick()*voices[j]:tick()
        end
        --outbuf[i] = outbuf[i] / total
    end
    ]]
    for i=0,frames-1 do        
        local t = lfo:tick()
        chorus:setModDepth(t)
        outbuf[i] = chorus:tick(clip:Tick(voicer:tick(),1.1,t))
    end
    
    for i=0,frames-1 do
        v[i+1] = outbuf[i]
    end
    diode:ProcessBlock(frames,v:data(),v:data())
    delay:ProcessBlock(frames,v:data(),v2:data())
    for i=0,frames-1 do
        outbuf[i] = v2[i+1] + v[i+1]
    end
end 

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end

function note_on(c,n,v)    
    voicer:noteOn(n,v)
end
--[[
function note_on(c,n,v)    
    local f = math.pow(2.0, (n-69)/12)*440.0            
    local found = false
    for i=1,8 do
        if(inuse[i] == false) then
            inuse[i] = true
            found = true
            envs[i]:keyOn()
            voices[i]:noteOn(f,v/127)
            freqs[i] = f
            vels[i] = v
            break
        end
    end
    if(found == false) then
        inuse[1] = true
        envs[1]:keyOn()
        voices[1]:noteOn(f,v/127)
        freqs[1] = f
        vels[1] = v
    end    
end
]]

function note_off(c,n,v)        
    voicer:noteOff(n,v)
end
--[[    
function note_off(c,n,v)        
    local f = math.pow(2.0, (n-69)/12)*440.0            
    for i=1,8 do
        if(freqs[i] == f) then 
            print(i)
            inuse[i] = false
            freqs[i] = 0
            vels[i] = 0
            envs[i]:keyOff()
            --voices[i]:noteOff(f)
            break
        end
    end
end
]]

function pitchbend(c,d1,d2)
    local x = 127*d2
    voicer:pitchBend(x)
end

function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 1) then
        voicer:controlChange(11,d2)
    elseif(d1 == 102) then        
        fc = d2/127                
    elseif(d1 == 103) then
        q = d2/127                
    end
end


function randomize()   
    for i=1,8 do     
        b3 = voices[i]
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
end
