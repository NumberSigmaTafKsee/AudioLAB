require('VCO')
require("VCF")
require("VCA")
require('audio_adsr')
require('stdsamples')
require('Amplifiers')
require('audio_simple_resampler')

adsr = audio_adsr.ADSR(0.1,0.2,0.8,0.2)
vco  = VCO.VCO(44100)
vcf  = VCF.VCF(44100.0,1000.0,0.5)
vca  = VCA.VCA(-6)
v    = stdsamples.sample_vector(256)
o    = stdsamples.sample_vector(512)
resmp= audio_simple_resampler.SimpleResampler()
resmp:setup(44100,2)

freq = 440
fc = 1.0
q  = 0.5
f0 = 261.6256
function freq2cv(f)   
    return math.log(f / f0);
end    
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
    vco:setFrequency(freq)    
    vco:ProcessSIMD(frames,v:data())        
    --[[
    resmp:up(frames,v:data(),o:data())    
    Amplifiers.udo1_simd(2*frames,o:data())
    Amplifiers.clamp_vector(2*frames,o:data(),o:data(),-1,1)
    resmp:down(frames,o:data(),v:data())        
    ]]
    vcf:setCutoff(fc*22050+fc*freq)
    vcf:setResonance(q)    
    vcf:ProcessSIMD(frames,v:data(),v:data())        
    vca:ProcessBlock(frames,v:data(),v:data())
    adsr:ProcessSIMD(frames,v:data(),v:data())
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
        fc = d2/127    
        fc = (math.pow(127,fc)-1.0)/126
    elseif(d1 == 103) then
        q = d2/127                
    end
end
