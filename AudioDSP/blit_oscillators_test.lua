require('va_blit_oscillators')
require('va_analog_svf')
require('audio_adsr')
require('stdsamples')
require('Amplifiers')

adsr = audio_adsr.ADSR(0.1,0.2,0.8,0.2)
osc  = va_blit_oscillators.blitSaw()
svf  = va_analog_svf.AnalogSVF(44100.0,1000.0,0.5)
v    = stdsamples.sample_vector(256)

freq = 0
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

function noise(input,output,frames)            
    local outbuf = new_buffer(output)             
    osc:ProcessSIMD(256,v:data())    
    Amplifiers.udo1_simd(256,v:data(),1.25)
    adsr:ProcessSIMD(256,v:data(),v:data())
    svf:setCutoff(fc)
    svf:setQ(q)
    svf:ProcessSIMD(256,v:data(),v:data())
    for i=0,frames-1 do
        outbuf[i] = v[i+1]
    end
end 

function freq_to_midi(f)
    return 12.0*math.log(f/440.0)/math.log(2) + 69
end 
function midi_to_freq(m)
    return math.pow(2.0, (m-69)/12)*440.0
end



function note_on(c,n,v)    
    freq = midi_to_freq(n)
    adsr:noteOn()
    osc:setFrequency(freq)
    
end

function note_off(c,n,v)        
    adsr:noteOff()
end

function pitchbend(c,d1,d2)
    local x = 127*d2
    
end

function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 1) then
        
    elseif(d1 == 102) then        
        fc = freq + 11025*d2/127                
    elseif(d1 == 103) then
        q = 20.0*d2/127                
    end
end
