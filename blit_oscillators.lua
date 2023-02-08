require('analog_blit_oscillators')
require('analog_krajeski_moog')
require('stdsamples')

saw = analog_blit_oscillators.blitSaw()
sqr = analog_blit_oscillators.blitSquare()
tri = analog_blit_oscillators.blitTriangle()
krajeski = analog_krajeski_moog.MoogLadder(44100)
krajeski:setDrive(2.0)

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
    local v = stdsamples.sample_vector(frames)
    saw:ProcessSIMD(frames,v:data())    
    krajeski:setFrequency(fc)
    krajeski:setQ(q)

    krajeski:ProcessSIMD(frames,v:data(),v:data())
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
    saw:setFrequency(midi_to_freq(n))
end

function note_off(c,n,v)        
    
end

function pitchbend(c,d1,d2)
    local x = 127*d2
end

function control(c,d1,d2)
    print(c,d1,d2)
    if(d1 == 1) then
    
    elseif(d1 == 102) then        
        fc = 22050/2*d2/127                
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
