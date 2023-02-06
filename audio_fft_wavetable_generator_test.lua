require('audiosystem')
require('audio_wavetable')
require('stk')

fft = require('audio_fft_wavetable_generator')
adsr = stk.ADSR()
adsr:setAttackRate(0)
adsr:setAllTimes(0.1,0.2,0.7,0.2)


freq = 440.0
numpress = 0
fc = 440.0
q  = 0.5


function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end


w = audio_wavetable.WaveTableOsc()
for i=1,129 do
    f = midi_to_freq(i)      
    v   = fft.FFTWaveTableGenerator.triangle(f,44100)    
    w:addWaveTable(v:size(),v:data(),f/44100)
end
w:setFrequency(440.0/44100)
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

function noise(input,output,frames)            
    local outbuf = new_buffer(output)            
    for i=0,frames-1 do        
        outbuf[i] = adsr:tick()*w:Tick()
    end
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
    w:setFrequency(freq/44100)    
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
